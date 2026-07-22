from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


def fix_one(sample_dir: Path) -> None:
    h5ad = sample_dir / "cell2loc.h5ad"
    if not h5ad.exists():
        return
    adata = ad.read_h5ad(h5ad)
    cell_types = list(adata.uns["mod"]["factor_names"])

    for key, csv_name in [("q05_cell_abundance_w_sf", "cell2loc_q05.csv"),
                          ("means_cell_abundance_w_sf", "cell2loc_mean.csv")]:
        arr = adata.obsm[key]
        if hasattr(arr, "values"):
            arr = arr.values
        df = pd.DataFrame(arr, index=adata.obs_names, columns=cell_types)
        df.to_csv(sample_dir / csv_name)
    print(f"  fixed: {sample_dir.name}")


def main() -> None:
    result_dir = Path("E:/GBM/results/c2l_anndata_ilc")
    dirs = sorted(d for d in result_dir.iterdir() if d.is_dir())
    for d in dirs:
        fix_one(d)
    print(f"done: {len(dirs)} samples")


if __name__ == "__main__":
    main()
