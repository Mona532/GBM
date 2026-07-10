"""Merge c2l results into Visium h5ad — output compatible with TLS pipeline"""
import argparse, json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

def merge_one(c2l_dir: Path, visium_dir: Path, out_path: Path) -> dict:
    sid = c2l_dir.name
    c2l_h5 = c2l_dir / "cell2loc.h5ad"
    if not c2l_h5.exists():
        return None

    # Load c2l posterior from already-run model output
    c2l = ad.read_h5ad(c2l_h5)
    q05_arr = c2l.obsm["q05_cell_abundance_w_sf"]
    mean_arr = c2l.obsm["means_cell_abundance_w_sf"]
    if hasattr(q05_arr, "values"):
        q05_arr = q05_arr.values
    if hasattr(mean_arr, "values"):
        mean_arr = mean_arr.values
    cell_types = list(c2l.uns["mod"]["factor_names"])

    # Load Visium raw data
    adata = sc.read_visium(visium_dir)
    adata.var_names_make_unique()
    # MT filter
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    adata = adata[:, keep].copy()

    # Align: c2l spot barcodes may differ from visium (in_tissue filtering)
    shared = adata.obs_names.intersection(c2l.obs_names)
    adata = adata[shared].copy()
    q05_arr = pd.DataFrame(q05_arr, index=c2l.obs_names, columns=cell_types).loc[shared].values
    mean_arr = pd.DataFrame(mean_arr, index=c2l.obs_names, columns=cell_types).loc[shared].values

    adata.obsm["c2l_ilc_q05"] = q05_arr
    adata.obsm["c2l_ilc_mean"] = mean_arr
    adata.uns["c2l_ilc_cell_types"] = np.array(cell_types, dtype=str)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_path, compression="gzip")
    print(f"  {sid}: {adata.n_obs} spots, {len(cell_types)} cell types")
    return {"sample_id": sid, "n_spots": adata.n_obs, "cell_types": len(cell_types)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c2l-dir", required=True, type=Path)
    parser.add_argument("--visium-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for d in sorted(args.c2l_dir.iterdir()):
        if not d.is_dir():
            continue
        sid = d.name
        visium = args.visium_dir / sid
        if not visium.exists():
            print(f"  SKIP {sid}: no visium dir")
            continue
        r = merge_one(d, visium, args.output / f"{sid}.h5ad")
        if r:
            rows.append(r)

    pd.DataFrame(rows).to_csv(args.output / "merge_summary.csv", index=False)
    print(f"Done: {len(rows)} samples")


if __name__ == "__main__":
    main()
