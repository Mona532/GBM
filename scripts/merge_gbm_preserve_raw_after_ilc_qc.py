from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix


ROOT = Path(__file__).resolve().parents[1]
QC_H5AD = ROOT / "ilc_qc_only" / "ilc_qc_only_filtered_raw.h5ad"
OUTPUT = ROOT / "gbm_merged_raw_counts_ilc_qc.h5ad"
SUMMARY = ROOT / "gbm_merged_raw_counts_ilc_qc_summary.md"


def find_first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError(f"None of these files exist: {paths}")


def read_local_10x(path: Path) -> ad.AnnData:
    matrix_path = find_first_existing([path / "matrix.mtx.gz", path / "matrix.mtx"])
    features_path = find_first_existing(
        [path / "features.tsv.gz", path / "features.tsv", path / "genes.tsv.gz", path / "genes.tsv"]
    )
    barcodes_path = find_first_existing([path / "barcodes.tsv.gz", path / "barcodes.tsv"])

    matrix = csr_matrix(mmread(matrix_path).T)
    features = pd.read_csv(features_path, sep="\t", header=None)
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)

    var = pd.DataFrame(index=features.iloc[:, 0].astype(str))
    var["gene_ids"] = features.iloc[:, 0].astype(str).values
    var["gene_symbols"] = features.iloc[:, 1].astype(str).values
    if features.shape[1] > 2:
        var["feature_types"] = features.iloc[:, 2].astype(str).values
    var.index.name = None

    obs = pd.DataFrame(index=barcodes.iloc[:, 0].astype(str))
    obs.index.name = None
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    if "barcode" not in meta.columns:
        raise ValueError(f"Missing barcode column in metadata: {path}")
    meta["barcode"] = meta["barcode"].astype(str)
    meta = meta.drop_duplicates(subset="barcode", keep="first")
    return meta.set_index("barcode", drop=False)


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    clean = df.copy()
    clean.index = pd.Index(clean.index.astype(str).tolist(), dtype=object)
    clean.index.name = None
    for column in clean.columns:
        series = clean[column]
        if pd.api.types.is_categorical_dtype(series):
            continue
        if pd.api.types.is_string_dtype(series) or pd.api.types.is_object_dtype(series):
            clean[column] = pd.Series(series.fillna("").astype(str).tolist(), index=clean.index, dtype=object)
    return clean


def load_ilc_raw(dataset: str, keep_barcodes: pd.Index) -> ad.AnnData:
    adata = read_local_10x(ROOT / dataset)
    adata.obs_names = adata.obs_names.astype(str)
    missing = keep_barcodes.difference(adata.obs_names)
    if len(missing) > 0:
        raise ValueError(f"{dataset} is missing {len(missing)} QC-selected barcodes")
    adata = adata[keep_barcodes].copy()
    adata.obs["barcode_raw"] = adata.obs_names
    adata.obs["dataset"] = dataset
    adata.obs["sample"] = dataset
    adata.obs["batch"] = dataset
    adata.obs_names = pd.Index([f"{dataset}_{x}" for x in adata.obs["barcode_raw"]], dtype=object)
    adata.obs["barcode_merged"] = adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    adata.obs = sanitize_dataframe(adata.obs)
    adata.var = sanitize_dataframe(adata.var)
    return adata


def load_rds_export(dataset: str) -> ad.AnnData:
    path = ROOT / f"{dataset}_10x"
    adata = read_local_10x(path)
    adata.obs_names = adata.obs_names.astype(str)
    meta = load_metadata(path / "metadata.csv")
    missing = adata.obs_names.difference(meta.index)
    if len(missing) > 0:
        raise ValueError(f"Metadata missing {len(missing)} barcodes for {dataset}")
    meta = meta.loc[adata.obs_names].copy()
    adata.obs["barcode_raw"] = adata.obs_names
    adata.obs["dataset"] = dataset
    adata.obs["sample"] = dataset
    adata.obs["batch"] = dataset
    adata.obs = adata.obs.join(meta.drop(columns=["barcode"], errors="ignore"))
    adata.obs_names = pd.Index([f"{dataset}_{x}" for x in adata.obs["barcode_raw"]], dtype=object)
    adata.obs["barcode_merged"] = adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    adata.obs = sanitize_dataframe(adata.obs)
    adata.var = sanitize_dataframe(adata.var)
    return adata


def write_summary(ilc_qc_obs: pd.DataFrame, merged: ad.AnnData) -> None:
    lines = [
        "# GBM merged raw-count object after ILC QC",
        "",
        "- Final matrix in `adata.X`: raw counts",
        "- Raw counts backup in `adata.layers['counts']`",
        "- `GBM01` and `GBM02`: unchanged",
        "- `GBM-ILC1` and `GBM-ILC2`: cells retained according to the QC skill output",
        "",
        "## ILC retention",
        "",
    ]
    for dataset in ["GBM-ILC1", "GBM-ILC2"]:
        subset = ilc_qc_obs.loc[ilc_qc_obs["sample"] == dataset]
        lines.append(f"- {dataset}: {subset.shape[0]} cells retained")
    lines.extend(["", "## Final object", "", f"- cells: {merged.n_obs}", f"- genes: {merged.n_vars}"])
    SUMMARY.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = True

    qc = ad.read_h5ad(QC_H5AD)
    qc_obs = qc.obs.copy()
    if "sample" not in qc_obs.columns or "barcode_raw" not in qc_obs.columns:
        raise ValueError("QC output is missing required obs columns: sample, barcode_raw")

    ilc1_keep = pd.Index(qc_obs.loc[qc_obs["sample"] == "GBM-ILC1", "barcode_raw"].astype(str))
    ilc2_keep = pd.Index(qc_obs.loc[qc_obs["sample"] == "GBM-ILC2", "barcode_raw"].astype(str))

    adatas = [
        load_ilc_raw("GBM-ILC1", ilc1_keep),
        load_ilc_raw("GBM-ILC2", ilc2_keep),
        load_rds_export("GBM01"),
        load_rds_export("GBM02"),
    ]
    merged = ad.concat(adatas, axis=0, join="outer", merge="first", index_unique=None, fill_value=0)
    merged.obs = sanitize_dataframe(merged.obs)
    merged.var = sanitize_dataframe(merged.var)
    merged.layers["counts"] = merged.X.copy()
    merged.uns["merge_note"] = (
        "GBM-ILC1 and GBM-ILC2 were QC filtered using the h5ad-qc-preannotation skill. "
        "GBM01 and GBM02 were kept unchanged. Final adata.X contains raw counts."
    )
    merged.write_h5ad(OUTPUT, compression="gzip")
    write_summary(qc_obs, merged)
    print(f"Wrote {OUTPUT}")
    print(f"Cells={merged.n_obs} Genes={merged.n_vars}")


if __name__ == "__main__":
    main()
