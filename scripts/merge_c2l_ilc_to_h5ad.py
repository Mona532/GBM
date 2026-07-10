from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def merge_one_sample(
    original_path: Path,
    q05_path: Path,
    mean_path: Path,
    output_path: Path,
) -> dict[str, object]:
    sample_id = original_path.stem
    print(f"[sample] {sample_id}")

    adata = ad.read_h5ad(original_path)
    q05 = pd.read_csv(q05_path, index_col=0)
    mean = pd.read_csv(mean_path, index_col=0)
    cell_types = list(q05.columns)

    # 只保留同时存在于 obsm 和 q05 中的 spot（in_tissue 筛选后可能不一致）
    shared_spots = adata.obs_names.intersection(q05.index)
    removed = len(q05) - len(shared_spots)
    if removed > 0:
        print(f"  removed {removed} spots not in original data")
    q05 = q05.loc[shared_spots]
    mean = mean.loc[shared_spots]

    adata.obsm["c2l_ilc_q05"] = q05.values
    adata.obsm["c2l_ilc_mean"] = mean.values
    adata.uns["c2l_ilc_cell_types"] = np.array(cell_types, dtype=str)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path, compression="gzip")
    print(f"  wrote: {adata.n_obs} spots, {len(cell_types)} cell types in obsm")
    return {"sample_id": sample_id, "n_spots": adata.n_obs, "cell_types": len(cell_types)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anndata-dir", required=True, type=Path)
    parser.add_argument("--c2l-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    result_dir = args.c2l_dir
    args.output.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample_dir in sorted(result_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_id = sample_dir.name
        original = args.anndata_dir / f"{sample_id}.h5ad"
        q05 = sample_dir / "cell2loc_q05.csv"
        mean = sample_dir / "cell2loc_mean.csv"

        if not original.exists():
            print(f"[skip] original h5ad missing: {sample_id}")
            continue
        if not q05.exists() or not mean.exists():
            print(f"[skip] c2l results missing: {sample_id}")
            continue

        rows.append(
            merge_one_sample(
                original_path=original,
                q05_path=q05,
                mean_path=mean,
                output_path=args.output / f"{sample_id}.h5ad",
            )
        )
        if args.limit and len(rows) >= args.limit:
            break

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output / "merge_summary.csv", index=False)
    print(f"[done] {len(rows)} samples")


if __name__ == "__main__":
    main()
