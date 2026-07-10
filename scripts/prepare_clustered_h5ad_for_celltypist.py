from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered.h5ad"
OUTPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered_symbolsafe_lognorm_for_celltypist.h5ad"


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
    if "lognorm" not in adata.layers:
        raise RuntimeError("Missing lognorm layer in clustered H5AD")
    adata.X = adata.layers["lognorm"].copy()
    adata = sanitize_for_h5ad_write(adata)
    adata.write_h5ad(OUTPUT, compression="gzip")
    print(f"written {OUTPUT}")
    print(f"shape {adata.shape}")


if __name__ == "__main__":
    main()
