"""Run cell2location on organized Visium spaceranger directories"""
from __future__ import annotations

import argparse, json
from pathlib import Path
import anndata as ad
import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import torch


def load_reference(ref_dir: Path) -> tuple[pd.DataFrame, pd.Index]:
    sig = pd.read_csv(ref_dir / "reference_signatures.csv", index_col=0)
    gene_txt = ref_dir / "selected_genes.txt"
    if gene_txt.exists():
        selected = pd.Index([x.strip() for x in gene_txt.read_text().splitlines() if x.strip()])
    else:
        selected = sig.index
    return sig, selected


def filter_mt_genes(adata: ad.AnnData) -> ad.AnnData:
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    return adata[:, keep].copy()


def run_one(sample_dir: Path, ref_sig: pd.DataFrame, selected: pd.Index,
            output_dir: Path, accel: str, max_epochs: int, batch_size: int) -> dict:
    sid = sample_dir.name
    sample_out = output_dir / sid
    out_h5ad = sample_out / "cell2loc.h5ad"
    out_csv = sample_out / "cell2loc_q05.csv"

    if out_h5ad.exists() and out_csv.exists():
        ck = pd.read_csv(out_csv, index_col=0)
        if ck.isna().mean().mean() < 0.01:
            return {"sample_id": sid, "n_spots": ck.shape[0], "n_cell_types": ck.shape[1], "status": "skipped"}

    sample_out.mkdir(parents=True, exist_ok=True)
    print(f"[sample] {sid}")

    adata = sc.read_visium(sample_dir)
    adata.var_names_make_unique()
    adata = filter_mt_genes(adata)

    shared = adata.var_names.intersection(ref_sig.index).intersection(selected)
    print(f"  genes shared: {len(shared)}")
    if len(shared) < 5000:
        print(f"  WARNING: very few shared genes!")
    adata = adata[:, shared].copy()
    sig = ref_sig.loc[shared, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata, batch_key=None)
    model = cell2location.models.Cell2location(
        adata, cell_state_df=sig, N_cells_per_location=8, detection_alpha=200,
    )
    model.train(max_epochs=max_epochs, accelerator=accel, batch_size=batch_size, train_size=1)
    model.save(sample_out / "model", overwrite=True)

    adata = model.export_posterior(adata, sample_kwargs={"num_samples": 1000, "batch_size": min(batch_size, adata.n_obs)})
    adata.write(out_h5ad, compression="gzip")

    cell_types = list(adata.uns["mod"]["factor_names"])
    for key, csv_name in [("q05_cell_abundance_w_sf", "cell2loc_q05.csv"),
                          ("means_cell_abundance_w_sf", "cell2loc_mean.csv")]:
        arr = adata.obsm[key]
        if hasattr(arr, "values"):
            arr = arr.values
        df = pd.DataFrame(arr, index=adata.obs_names, columns=cell_types)
        df.to_csv(sample_out / csv_name)

    with open(sample_out / "meta.json", "w") as f:
        json.dump({"sample_id": sid, "n_spots": int(adata.n_obs), "n_genes": int(adata.n_vars),
                   "cell_types": cell_types, "max_epochs": max_epochs, "batch_size": batch_size}, f, indent=2)
    print(f"  done: {adata.n_obs} spots, {len(cell_types)} cell types")
    return {"sample_id": sid, "n_spots": int(adata.n_obs), "n_cell_types": len(cell_types), "status": "done"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--visium-dir", required=True, type=Path)
    parser.add_argument("--ref-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-epochs", type=int, default=2000)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    torch.set_float32_matmul_precision("medium")
    accel = "gpu" if torch.cuda.is_available() else "cpu"
    args.output.mkdir(parents=True, exist_ok=True)
    ref_sig, selected = load_reference(args.ref_dir)

    dirs = sorted(p for p in args.visium_dir.iterdir() if p.is_dir() and (p / "filtered_feature_bc_matrix.h5").exists())
    if args.limit > 0:
        dirs = dirs[:args.limit]

    print(f"Samples: {len(dirs)} | Ref types: {ref_sig.shape[1]} | Genes: {len(selected)} | Accel: {accel}")

    rows = []
    for d in dirs:
        rows.append(run_one(d, ref_sig, selected, args.output, accel, args.max_epochs, args.batch_size))

    pd.DataFrame(rows).to_csv(args.output / "summary.csv", index=False)
    done = sum(1 for r in rows if r["status"] == "done")
    skip = sum(1 for r in rows if r["status"] == "skipped")
    print(f"[done] {done} processed, {skip} skipped")


if __name__ == "__main__":
    main()
