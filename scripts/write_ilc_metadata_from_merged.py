from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_merged_raw_counts_ilc_qc.h5ad"


def write_dataset_metadata(adata, dataset: str, out_csv: Path) -> None:
    sub = adata.obs.loc[adata.obs["dataset"].astype(str) == dataset].copy()
    sub["barcode"] = sub["barcode_raw"].astype(str)
    cols = ["barcode"] + [c for c in sub.columns if c != "barcode"]
    sub = sub.loc[:, cols]
    sub.to_csv(out_csv, index=False)


def main() -> None:
    adata = sc.read_h5ad(INPUT)
    for dataset in ["GBM-ILC1", "GBM-ILC2"]:
        out_csv = ROOT / dataset / "metadata.csv"
        write_dataset_metadata(adata, dataset, out_csv)
        print(f"written {out_csv}")


if __name__ == "__main__":
    main()
