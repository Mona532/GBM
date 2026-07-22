"""Merge consolidated C2L results into h5ad obsm"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

C2L = Path(r"E:/GBM/results/c2l_consolidated")
H5AD = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata")
OUT = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated")
OUT.mkdir(parents=True, exist_ok=True)

for c2l_dir in sorted(C2L.iterdir()):
    if not c2l_dir.is_dir(): continue
    sid = c2l_dir.name
    q05_csv = c2l_dir / "cell2loc_q05.csv"
    mean_csv = c2l_dir / "cell2loc_mean.csv"
    h5 = H5AD / f"{sid}.h5ad"
    if not q05_csv.exists() or not h5.exists(): continue

    q05 = pd.read_csv(q05_csv, index_col=0)
    mean = pd.read_csv(mean_csv, index_col=0)
    adata = ad.read_h5ad(h5)
    shared = adata.obs_names.intersection(q05.index)
    if len(shared) < 100: continue
    adata = adata[shared]
    q05 = q05.loc[shared]; mean = mean.loc[shared]
    adata.obsm["c2l_ilc_q05"] = q05.values
    adata.obsm["c2l_ilc_mean"] = mean.values
    adata.uns["c2l_ilc_cell_types"] = np.array(q05.columns, dtype=str)
    adata.write(OUT / f"{sid}.h5ad", compression="gzip")
    if int(sid[-1]) == 1 and len(list(C2L.iterdir())) > 90:  # print only first few
        print(f"  {sid}: {adata.n_obs} spots, {q05.shape[1]} types")

print(f"Done: {len(list(OUT.glob('*.h5ad')))} h5ad files")
