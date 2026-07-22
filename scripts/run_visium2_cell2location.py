from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch


def set_var_names_from_symbol(adata: ad.AnnData, symbol_col: str = "SYMBOL") -> ad.AnnData:
    if symbol_col in adata.var.columns:
        symbols = adata.var[symbol_col].astype(str).to_numpy()
        fallback = adata.var_names.astype(str).to_numpy()
        use_names = np.where((symbols != "nan") & (symbols != ""), symbols, fallback)
        adata.var_names = pd.Index(use_names)
    adata.var_names_make_unique()
    return adata


def ensure_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    return adata


def filter_mt_genes(adata: ad.AnnData) -> ad.AnnData:
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    return adata[:, keep].copy()


def subset_gene_expression_only(adata: ad.AnnData) -> ad.AnnData:
    if "feature_types" in adata.var.columns:
        keep = adata.var["feature_types"].astype(str).eq("Gene Expression").to_numpy()
        adata = adata[:, keep].copy()
    return adata


def balanced_subsample(
    adata: ad.AnnData,
    label_key: str,
    max_cells_per_label: int,
    seed: int,
) -> ad.AnnData:
    if max_cells_per_label <= 0:
        return adata
    keep_idx: list[np.ndarray] = []
    rng = np.random.default_rng(seed)
    labels = adata.obs[label_key].astype(str)
    for label, idx in labels.groupby(labels).groups.items():
        idx = np.asarray(list(idx))
        if len(idx) > max_cells_per_label:
            idx = rng.choice(idx, size=max_cells_per_label, replace=False)
        keep_idx.append(idx)
    keep_idx = np.concatenate(keep_idx)
    return adata[keep_idx].copy()


def train_reference_model(
    ref_path: Path,
    out_dir: Path,
    label_key: str,
    batch_key: str | None,
    exclude_labels: list[str],
    max_cells_per_label: int,
    seed: int,
    accelerator: str,
) -> tuple[ad.AnnData, pd.DataFrame]:
    print("=" * 80)
    print("Loading snRNA reference")
    print("=" * 80)
    posterior_path = out_dir / "reference_posterior.h5ad"
    signatures_path = out_dir / "reference_signatures.csv"
    model_dir = out_dir / "reference_model"
    if posterior_path.exists() and model_dir.exists():
        print("Reusing existing reference model and posterior")
        adata_ref = ad.read_h5ad(posterior_path)
        if "means_per_cluster_mu_fg" not in adata_ref.varm:
            raise RuntimeError("reference_posterior.h5ad exists but lacks means_per_cluster_mu_fg")
        inf_aver = extract_inf_aver(adata_ref)
        inf_aver.to_csv(signatures_path)
        return adata_ref, inf_aver

    adata_ref = ad.read_h5ad(ref_path)
    adata_ref = ensure_raw_counts(adata_ref)
    adata_ref = set_var_names_from_symbol(adata_ref, symbol_col="SYMBOL")
    adata_ref = filter_mt_genes(adata_ref)

    if label_key not in adata_ref.obs.columns:
        raise KeyError(f"Reference label_key '{label_key}' not found in obs")

    adata_ref = adata_ref[~adata_ref.obs[label_key].astype(str).isin(exclude_labels)].copy()
    sc.pp.filter_genes(adata_ref, min_cells=3)
    adata_ref = balanced_subsample(adata_ref, label_key, max_cells_per_label, seed)

    print(f"Reference after filtering: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")
    print(adata_ref.obs[label_key].astype(str).value_counts().to_string())

    kwargs = {"labels_key": label_key}
    if batch_key and batch_key in adata_ref.obs.columns:
        kwargs["batch_key"] = batch_key

    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, **kwargs)
    mod = cell2location.models.RegressionModel(adata_ref)
    mod.train(max_epochs=150, batch_size=4096, accelerator=accelerator)

    ref_model_dir = out_dir / "reference_model"
    mod.save(ref_model_dir, overwrite=True)

    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": 1000, "batch_size": 4096},
    )
    adata_ref.write(out_dir / "reference_posterior.h5ad")

    inf_aver = extract_inf_aver(adata_ref)
    inf_aver.to_csv(out_dir / "reference_signatures.csv")
    return adata_ref, inf_aver


def extract_inf_aver(adata_ref: ad.AnnData) -> pd.DataFrame:
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"]
    if isinstance(inf_aver, pd.DataFrame):
        inf_aver = inf_aver.copy()
    else:
        inf_aver = pd.DataFrame(inf_aver, index=adata_ref.var_names)
    inf_aver.index = adata_ref.var_names
    inf_aver.columns = [str(c).replace("means_per_cluster_mu_fg_", "") for c in inf_aver.columns]
    return inf_aver


def run_spatial_model(
    spatial_path: Path,
    out_dir: Path,
    inf_aver: pd.DataFrame,
    accelerator: str,
    max_epochs: int,
    batch_size: int,
) -> ad.AnnData:
    print("=" * 80)
    print("Loading integrated Visium object")
    print("=" * 80)
    adata_vis = ad.read_h5ad(spatial_path)
    adata_vis = subset_gene_expression_only(adata_vis)
    adata_vis = filter_mt_genes(adata_vis)
    adata_vis.var_names_make_unique()

    shared = adata_vis.var_names.intersection(inf_aver.index)
    if len(shared) < 1000:
        raise RuntimeError(f"Too few shared genes between reference and spatial data: {len(shared)}")
    print(f"Spatial object: {adata_vis.n_obs} spots x {adata_vis.n_vars} genes")
    print(f"Shared genes: {len(shared)}")

    adata_vis = adata_vis[:, shared].copy()
    inf_aver = inf_aver.loc[shared, :].copy()

    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis,
        batch_key="sample_id" if "sample_id" in adata_vis.obs.columns else None,
    )
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=8,
        detection_alpha=20,
    )
    mod.train(
        max_epochs=max_epochs,
        accelerator=accelerator,
        batch_size=batch_size,
        train_size=1,
    )
    spatial_model_dir = out_dir / "spatial_model"
    mod.save(spatial_model_dir, overwrite=True)

    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={"num_samples": 1000, "batch_size": adata_vis.n_obs},
    )
    adata_vis.write(out_dir / "visium2_cell2location.h5ad", compression="gzip")

    cell_types = list(adata_vis.uns["mod"]["factor_names"])
    q05 = pd.DataFrame(
        adata_vis.obsm["q05_cell_abundance_w_sf"],
        index=adata_vis.obs_names,
        columns=cell_types,
    )
    means = pd.DataFrame(
        adata_vis.obsm["means_cell_abundance_w_sf"],
        index=adata_vis.obs_names,
        columns=cell_types,
    )
    q05.to_csv(out_dir / "visium2_cell2location_q05.csv")
    means.to_csv(out_dir / "visium2_cell2location_mean.csv")

    obs_export = adata_vis.obs.copy()
    for col in q05.columns:
        obs_export[f"cell2loc_q05__{col}"] = q05[col]
    obs_export.to_csv(out_dir / "visium2_obs_with_cell2loc_q05.csv")
    return adata_vis


def main() -> None:
    parser = argparse.ArgumentParser(description="Run cell2location on integrated Visium 2 dataset.")
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--spatial", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--label-key", default="annotation_coarse")
    parser.add_argument("--batch-key", default="donor_id")
    parser.add_argument("--exclude-label", action="append", default=["Undefined"])
    parser.add_argument("--max-cells-per-label", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--spatial-max-epochs", type=int, default=20000)
    parser.add_argument("--spatial-batch-size", type=int, default=2048)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_float32_matmul_precision("medium")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"Using accelerator: {accelerator}")

    _, inf_aver = train_reference_model(
        ref_path=args.reference,
        out_dir=args.output,
        label_key=args.label_key,
        batch_key=args.batch_key,
        exclude_labels=args.exclude_label,
        max_cells_per_label=args.max_cells_per_label,
        seed=args.seed,
        accelerator=accelerator,
    )
    adata_vis = run_spatial_model(
        spatial_path=args.spatial,
        out_dir=args.output,
        inf_aver=inf_aver,
        accelerator=accelerator,
        max_epochs=args.spatial_max_epochs,
        batch_size=args.spatial_batch_size,
    )

    meta = {
        "reference": str(args.reference),
        "spatial": str(args.spatial),
        "label_key": args.label_key,
        "batch_key": args.batch_key,
        "exclude_labels": args.exclude_label,
        "max_cells_per_label": args.max_cells_per_label,
        "seed": args.seed,
        "accelerator": accelerator,
        "spatial_max_epochs": args.spatial_max_epochs,
        "spatial_batch_size": args.spatial_batch_size,
        "n_spots": int(adata_vis.n_obs),
        "n_genes": int(adata_vis.n_vars),
        "cell_types": list(adata_vis.uns["mod"]["factor_names"]),
    }
    with open(args.output / "run_config.json", "w", encoding="utf-8") as handle:
        json.dump(meta, handle, ensure_ascii=False, indent=2)

    print("[done]")


if __name__ == "__main__":
    main()
