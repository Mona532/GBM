from __future__ import annotations

import gzip
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_merged_raw_counts_ilc_qc.h5ad"
GTF = Path(r"E:\blood_cell_development\geo_data\geo_bundle\gtf\Homo_sapiens.GRCh38.115.gtf.gz")
OUTPUT = ROOT / "gbm_merged_raw_counts_ilc_qc_symbol_safe.h5ad"
SUMMARY = ROOT / "gbm_merged_raw_counts_ilc_qc_symbol_safe_summary.md"


def load_gtf_map(gtf_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attrs)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attrs)
            if not gene_id_match or not gene_name_match:
                continue
            gene_id = gene_id_match.group(1).split(".")[0]
            gene_name = gene_name_match.group(1)
            if gene_id and gene_name and gene_id not in mapping:
                mapping[gene_id] = gene_name
    return mapping


def aggregate_duplicate_columns(matrix: sp.csr_matrix, names: list[str]) -> tuple[sp.csr_matrix, list[str]]:
    name_to_idx: dict[str, int] = {}
    col_idx = []
    unique_names = []
    for name in names:
        if name not in name_to_idx:
            name_to_idx[name] = len(unique_names)
            unique_names.append(name)
        col_idx.append(name_to_idx[name])
    gene_idx = np.arange(len(names), dtype=int)
    data = np.ones(len(names), dtype=np.float32)
    indicator = sp.csr_matrix((data, (gene_idx, np.array(col_idx))), shape=(len(names), len(unique_names)))
    return (matrix @ indicator).tocsr(), unique_names


def sanitize_for_h5ad_write(adata: ad.AnnData) -> ad.AnnData:
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = True
    adata.obs_names = pd.Index(adata.obs_names.astype(str))
    adata.var_names = pd.Index(adata.var_names.astype(str))
    for frame in [adata.obs, adata.var]:
        frame.index = pd.Index(frame.index.astype(str))
        for col in frame.columns:
            if isinstance(frame[col].dtype, pd.CategoricalDtype):
                continue
            if pd.api.types.is_string_dtype(frame[col]):
                frame[col] = frame[col].fillna("").astype(str).astype(object)
    return adata


def main() -> None:
    adata = sc.read_h5ad(INPUT)
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    gene_map = load_gtf_map(GTF)

    original_var_names = adata.var_names.astype(str)
    original_gene_ids = adata.var["gene_ids"].astype(str) if "gene_ids" in adata.var.columns else original_var_names

    target_symbols: list[str] = []
    source_types: list[str] = []
    mapped_from_ensembl = 0
    unmapped_ensembl = 0
    for var_name, gene_id in zip(original_var_names, original_gene_ids):
        gene_id_base = gene_id.split(".")[0]
        if re.match(r"^ENSG\d+", gene_id_base):
            symbol = gene_map.get(gene_id_base)
            if symbol:
                target_symbols.append(symbol)
                source_types.append("ensembl_to_symbol")
                mapped_from_ensembl += 1
            else:
                target_symbols.append(gene_id_base)
                source_types.append("ensembl_unmapped")
                unmapped_ensembl += 1
        else:
            target_symbols.append(var_name)
            source_types.append("symbol_passthrough")

    duplicate_before = int(pd.Index(target_symbols).duplicated().sum())
    new_counts, unique_symbols = aggregate_duplicate_columns(counts, target_symbols)

    var_meta = pd.DataFrame(
        {
            "target_symbol": target_symbols,
            "original_var_name": original_var_names,
            "original_gene_id": original_gene_ids,
            "symbol_source_type": source_types,
        }
    )
    collapsed_meta = (
        var_meta.groupby("target_symbol", sort=False)
        .agg(
            original_var_name=("original_var_name", lambda x: "|".join(pd.unique(pd.Series(x).astype(str))[:10])),
            original_gene_id=("original_gene_id", lambda x: "|".join(pd.unique(pd.Series(x).astype(str))[:10])),
            symbol_source_type=("symbol_source_type", lambda x: "|".join(pd.unique(pd.Series(x).astype(str)))),
        )
        .reindex(unique_symbols)
    )
    collapsed_meta.index = pd.Index(unique_symbols)
    collapsed_meta["gene_symbols"] = collapsed_meta.index.astype(str)
    collapsed_meta["gene_ids"] = collapsed_meta["original_gene_id"].astype(str)

    new = ad.AnnData(X=new_counts, obs=adata.obs.copy(), var=collapsed_meta)
    new.var_names = pd.Index(unique_symbols)
    new.layers["counts"] = new_counts.copy()
    if "lognorm" in adata.layers:
        new.layers["lognorm"] = None  # incompatible after collapse; do not carry stale layer
        del new.layers["lognorm"]
    for key in ["uns", "obsm", "obsp", "varm"]:
        if hasattr(adata, key):
            getattr(new, key).update(getattr(adata, key))

    summary = [
        "# Gene Standardization Summary",
        "",
        f"- input genes: {adata.n_vars}",
        f"- output genes: {new.n_vars}",
        f"- ensembl mapped to symbol: {mapped_from_ensembl}",
        f"- ensembl unmapped: {unmapped_ensembl}",
        f"- duplicate symbols collapsed: {duplicate_before}",
        f"- gtf: `{GTF}`",
    ]
    SUMMARY.write_text("\n".join(summary), encoding="utf-8")

    new = sanitize_for_h5ad_write(new)
    new.write_h5ad(OUTPUT, compression="gzip")
    print(f"written {OUTPUT}")
    print(f"shape {new.shape}")


if __name__ == "__main__":
    main()
