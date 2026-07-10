from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import torch


def ensure_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    return adata.raw.to_adata() if adata.raw is not None else adata


def filter_mt_genes(adata: ad.AnnData) -> ad.AnnData:
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    return adata[:, keep].copy()


def extract_inf_aver(adata_ref: ad.AnnData) -> pd.DataFrame:
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"]
    if isinstance(inf_aver, pd.DataFrame):
        inf_aver = inf_aver.copy()
    else:
        inf_aver = pd.DataFrame(inf_aver, index=adata_ref.var_names)
    inf_aver.index = adata_ref.var_names
    inf_aver.columns = [str(c).replace("means_per_cluster_mu_fg_", "") for c in inf_aver.columns]
    return inf_aver


def main() -> None:
    parser = argparse.ArgumentParser(description="Train cell2location RegressionModel and export signatures only.")
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--label-key", required=True)
    parser.add_argument("--batch-key", default="")
    parser.add_argument("--max-epochs", type=int, default=150)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--posterior-batch-size", type=int, default=4096)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_float32_matmul_precision("medium")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"

    adata_ref = ad.read_h5ad(args.reference)
    adata_ref = ensure_raw_counts(adata_ref)
    adata_ref.var_names_make_unique()
    adata_ref = filter_mt_genes(adata_ref)
    sc.pp.filter_genes(adata_ref, min_cells=3)

    kwargs = {"labels_key": args.label_key}
    if args.batch_key and args.batch_key in adata_ref.obs.columns:
        kwargs["batch_key"] = args.batch_key

    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, **kwargs)
    mod = cell2location.models.RegressionModel(adata_ref)
    mod.train(max_epochs=args.max_epochs, batch_size=args.batch_size, accelerator=accelerator)

    model_dir = args.output / "reference_model"
    mod.save(model_dir, overwrite=True)
    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": 1000, "batch_size": args.posterior_batch_size},
    )
    posterior_path = args.output / "reference_posterior.h5ad"
    signatures_path = args.output / "reference_signatures.csv"
    adata_ref.write(posterior_path)
    inf_aver = extract_inf_aver(adata_ref)
    inf_aver.to_csv(signatures_path)

    with open(args.output / "reference_train_meta.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "reference": str(args.reference),
                "label_key": args.label_key,
                "batch_key": args.batch_key,
                "accelerator": accelerator,
                "max_epochs": args.max_epochs,
                "batch_size": args.batch_size,
                "posterior_batch_size": args.posterior_batch_size,
                "n_cells": int(adata_ref.n_obs),
                "n_genes": int(adata_ref.n_vars),
                "cell_types": list(inf_aver.columns),
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )
    print("[done]")


if __name__ == "__main__":
    main()
