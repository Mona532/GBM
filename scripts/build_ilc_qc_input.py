from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "ilc_qc_input_raw.h5ad"

DATASETS = [
    ("GBM-ILC1", ROOT / "GBM-ILC1"),
    ("GBM-ILC2", ROOT / "GBM-ILC2"),
]


def load_10x(dataset: str, path: Path) -> ad.AnnData:
    adata = sc.read_10x_mtx(path, var_names="gene_ids", make_unique=True, gex_only=True)
    adata.obs_names = adata.obs_names.astype(str)
    adata.obs["barcode_raw"] = adata.obs_names
    adata.obs["dataset"] = dataset
    adata.obs["sample"] = dataset
    adata.obs["batch"] = dataset
    adata.obs_names = pd.Index([f"{dataset}_{x}" for x in adata.obs["barcode_raw"]], dtype=object)
    adata.obs["barcode_merged"] = adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    return adata


def main() -> None:
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = True
    adatas = [load_10x(dataset, path) for dataset, path in DATASETS]
    merged = ad.concat(adatas, axis=0, join="outer", merge="same", index_unique=None, fill_value=0)
    merged.layers["counts"] = merged.X.copy()
    merged.write_h5ad(OUTPUT, compression="gzip")
    print(f"Wrote {OUTPUT}")
    print(f"Cells={merged.n_obs} Genes={merged.n_vars}")


if __name__ == "__main__":
    main()
