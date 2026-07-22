from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scrublet as scr
from scipy.io import mmread
from scipy.sparse import csr_matrix
from anndata import AnnData, concat
from pandas.api.types import (
    is_bool_dtype,
    is_categorical_dtype,
    is_integer_dtype,
    is_float_dtype,
    is_object_dtype,
    is_string_dtype,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "gbm_merged_ilc_qc_counts_metadata.h5ad"
SUMMARY = ROOT / "gbm_merged_ilc_qc_summary.md"

DATASETS = [
    {
        "name": "GBM-ILC1",
        "path": ROOT / "GBM-ILC1",
        "metadata": None,
        "sample": "GBM-ILC1",
        "apply_qc": True,
    },
    {
        "name": "GBM-ILC2",
        "path": ROOT / "GBM-ILC2",
        "metadata": None,
        "sample": "GBM-ILC2",
        "apply_qc": True,
    },
    {
        "name": "GBM01",
        "path": ROOT / "GBM01_10x",
        "metadata": ROOT / "GBM01_10x" / "metadata.csv",
        "sample": "GBM01",
        "apply_qc": False,
    },
    {
        "name": "GBM02",
        "path": ROOT / "GBM02_10x",
        "metadata": ROOT / "GBM02_10x" / "metadata.csv",
        "sample": "GBM02",
        "apply_qc": False,
    },
]


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    if "barcode" not in meta.columns:
        raise ValueError(f"Missing barcode column in metadata: {path}")
    meta["barcode"] = meta["barcode"].astype(str)
    meta = meta.drop_duplicates(subset="barcode", keep="first")
    return meta.set_index("barcode", drop=False)


def find_first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError(f"None of these files exist: {paths}")


def read_local_10x(path: Path) -> AnnData:
    matrix_path = find_first_existing([path / "matrix.mtx.gz", path / "matrix.mtx"])
    features_path = find_first_existing(
        [
            path / "features.tsv.gz",
            path / "features.tsv",
            path / "genes.tsv.gz",
            path / "genes.tsv",
        ]
    )
    barcodes_path = find_first_existing([path / "barcodes.tsv.gz", path / "barcodes.tsv"])

    matrix = csr_matrix(mmread(matrix_path).T)
    features = pd.read_csv(features_path, sep="\t", header=None)
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)

    if features.shape[1] < 2:
        raise ValueError(f"Unexpected features format: {features_path}")

    var = pd.DataFrame(index=features.iloc[:, 0].astype(str))
    var["gene_ids"] = features.iloc[:, 0].astype(str).values
    var["gene_symbols"] = features.iloc[:, 1].astype(str).values
    if features.shape[1] > 2:
        var["feature_types"] = features.iloc[:, 2].astype(str).values
    var.index.name = None

    obs = pd.DataFrame(index=barcodes.iloc[:, 0].astype(str))
    obs.index.name = None
    adata = AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    clean = df.copy()
    clean.index = pd.Index(clean.index.astype(str).tolist(), dtype=object)
    clean.index.name = None
    for column in clean.columns:
        series = clean[column]
        if is_categorical_dtype(series) or is_bool_dtype(series) or is_integer_dtype(series) or is_float_dtype(series):
            continue
        if is_string_dtype(series) or is_object_dtype(series):
            clean[column] = pd.Series(series.astype(str).tolist(), index=clean.index, dtype=object)
    return clean


def mad(x: pd.Series | np.ndarray) -> float:
    arr = np.asarray(x, dtype=float)
    median = np.median(arr)
    return float(np.median(np.abs(arr - median)))


def qc_filter_10x(adata: AnnData, dataset_name: str) -> tuple[AnnData, dict]:
    qc = adata.copy()
    gene_symbols = qc.var["gene_symbols"].copy()
    gene_symbols = gene_symbols.where(gene_symbols.notna(), pd.Series(qc.var_names, index=qc.var_names))
    qc.var["mt"] = gene_symbols.astype(str).str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(qc, qc_vars=["mt"], inplace=True)

    log_genes = np.log1p(qc.obs["n_genes_by_counts"].to_numpy())
    log_counts = np.log1p(qc.obs["total_counts"].to_numpy())
    mt_pct = qc.obs["pct_counts_mt"].to_numpy()

    gene_floor = max(np.expm1(np.median(log_genes) - 3 * mad(log_genes)), 200.0)
    count_floor = max(np.expm1(np.median(log_counts) - 3 * mad(log_counts)), 500.0)
    mt_ceiling = min(np.median(mt_pct) + 3 * mad(mt_pct), 40.0)

    low_quality_mask = (
        (qc.obs["n_genes_by_counts"] < gene_floor)
        | (qc.obs["total_counts"] < count_floor)
        | (qc.obs["pct_counts_mt"] > mt_ceiling)
    )

    qc.obs["low_quality"] = low_quality_mask.to_numpy(dtype=bool)
    passing = qc[~qc.obs["low_quality"]].copy()

    scrub = scr.Scrublet(passing.X, expected_doublet_rate=0.06, random_state=0)
    scores, predicted = scrub.scrub_doublets(
        use_approx_neighbors=True,
        min_counts=3,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30,
        verbose=False,
    )
    passing.obs["doublet_score"] = scores
    passing.obs["predicted_doublet"] = predicted.astype(bool)

    qc.obs["doublet_score"] = np.nan
    qc.obs["predicted_doublet"] = False
    qc.obs.loc[passing.obs_names, "doublet_score"] = passing.obs["doublet_score"].to_numpy()
    qc.obs.loc[passing.obs_names, "predicted_doublet"] = passing.obs["predicted_doublet"].to_numpy()

    doublet_mask = qc.obs["predicted_doublet"].fillna(False).to_numpy(dtype=bool)
    keep_mask = ~(qc.obs["low_quality"].to_numpy(dtype=bool) | doublet_mask)

    filtered = qc[keep_mask].copy()
    summary = {
        "dataset": dataset_name,
        "input_cells": int(qc.n_obs),
        "kept_cells": int(filtered.n_obs),
        "removed_low_quality": int(qc.obs["low_quality"].sum()),
        "removed_doublets": int(doublet_mask.sum()),
        "gene_floor": float(gene_floor),
        "count_floor": float(count_floor),
        "mt_ceiling": float(mt_ceiling),
    }
    return filtered, summary


def load_dataset(spec: dict) -> AnnData:
    adata = read_local_10x(spec["path"])

    adata.obs_names = adata.obs_names.astype(str)
    adata.obs["barcode_raw"] = adata.obs_names
    adata.obs["dataset"] = spec["name"]
    adata.obs["sample"] = spec["sample"]
    adata.obs["batch"] = spec["name"]

    if spec["metadata"] is not None:
        meta = load_metadata(spec["metadata"])
        missing = adata.obs_names.difference(meta.index)
        if len(missing) > 0:
            raise ValueError(
                f"Metadata missing {len(missing)} barcodes for {spec['name']}"
            )
        meta = meta.loc[adata.obs_names].copy()
        adata.obs = adata.obs.join(meta.drop(columns=["barcode"], errors="ignore"))

    adata.obs_names = [f"{spec['name']}_{barcode}" for barcode in adata.obs["barcode_raw"]]
    adata.obs["barcode_merged"] = adata.obs_names

    qc_summary = None
    if spec.get("apply_qc", False):
        adata, qc_summary = qc_filter_10x(adata, spec["name"])
        adata.obs["qc_applied"] = True
    else:
        adata.obs["qc_applied"] = False

    adata.obs = sanitize_dataframe(adata.obs)
    adata.var = sanitize_dataframe(adata.var)
    return adata, qc_summary


def write_summary(qc_summaries: list[dict], merged: AnnData) -> None:
    lines = [
        "# GBM merge summary",
        "",
        "- Output: `gbm_merged_ilc_qc_counts_metadata.h5ad`",
        "- Merge rule: gene union, fill missing genes with 0",
        "- Raw counts retained in `adata.X` and `adata.layers['counts']`",
        "- `GBM01` and `GBM02`: no QC applied",
        "- `GBM-ILC1` and `GBM-ILC2`: low-quality cells and Scrublet-predicted doublets removed",
        "",
        "## Dataset QC",
        "",
    ]
    for item in qc_summaries:
        if item is None:
            continue
        lines.extend(
            [
                f"### {item['dataset']}",
                "",
                f"- input_cells: {item['input_cells']}",
                f"- kept_cells: {item['kept_cells']}",
                f"- removed_low_quality: {item['removed_low_quality']}",
                f"- removed_doublets: {item['removed_doublets']}",
                f"- gene_floor: {item['gene_floor']:.2f}",
                f"- count_floor: {item['count_floor']:.2f}",
                f"- mt_ceiling: {item['mt_ceiling']:.2f}",
                "",
            ]
        )
    lines.extend(
        [
            "## Final object",
            "",
            f"- cells: {merged.n_obs}",
            f"- genes: {merged.n_vars}",
        ]
    )
    SUMMARY.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ad.settings.allow_write_nullable_strings = True
    loaded = [load_dataset(spec) for spec in DATASETS]
    adatas = [item[0] for item in loaded]
    qc_summaries = [item[1] for item in loaded]
    merged = concat(
        adatas,
        axis=0,
        join="outer",
        merge="first",
        label=None,
        index_unique=None,
        fill_value=0,
    )
    merged.obs = sanitize_dataframe(merged.obs)
    merged.var = sanitize_dataframe(merged.var)
    merged.layers["counts"] = merged.X.copy()
    merged.uns["merge_note"] = (
        "Merged 4 GBM datasets by gene union. "
        "GBM-ILC1 and GBM-ILC2 received low-quality-cell filtering and Scrublet doublet removal. "
        "GBM01 and GBM02 were not filtered. "
        "adata.X and adata.layers['counts'] both contain raw counts."
    )
    merged.write_h5ad(OUTPUT)
    write_summary(qc_summaries, merged)
    print(f"Wrote {OUTPUT}")
    print(f"Cells={merged.n_obs} Genes={merged.n_vars}")


if __name__ == "__main__":
    main()
