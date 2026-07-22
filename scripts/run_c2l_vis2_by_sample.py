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
from cell2location.utils.filtering import filter_genes


REGION_MAP = {
    "C": "Cortex",
    "T": "Tumor",
    "TC": "TumorCore",
    "TI": "TumorInfiltration",
}


def parse_sample_metadata(sample_name: str) -> dict[str, str]:
    clean = sample_name.lstrip("#")
    parts = clean.split("_")
    region_code = parts[-2]
    return {
        "sample_id": clean,
        "sample_dir_name": sample_name,
        "patient_id": parts[0],
        "region_code": region_code,
        "region_label": REGION_MAP.get(region_code, region_code),
        "idh_status": "IDHmutant" if "IDHMutant" in clean else "IDHunknown",
        "cohort": "UniversityClinicFreiburg",
    }


def filter_mt_genes(adata: ad.AnnData) -> ad.AnnData:
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    return adata[:, keep].copy()


def load_reference_signatures(ref_dir: Path) -> pd.DataFrame:
    sig = pd.read_csv(ref_dir / "reference_signatures.csv", index_col=0)
    if (sig.sum(axis=0) <= 0).any():
        raise RuntimeError("reference_signatures.csv contains zero-sum columns")
    return sig


def load_selected_genes(ref_dir: Path) -> pd.Index:
    genes_txt = ref_dir / "selected_genes.txt"
    if genes_txt.exists():
        genes = [x.strip() for x in genes_txt.read_text(encoding="utf-8").splitlines() if x.strip()]
        return pd.Index(genes)

    ref_posterior = ad.read_h5ad(ref_dir / "reference_posterior.h5ad")
    genes = filter_genes(
        ref_posterior,
        cell_count_cutoff=15,
        cell_percentage_cutoff2=0.05,
        nonz_mean_cutoff=1.12,
    )
    genes_txt.write_text("\n".join(map(str, genes)), encoding="utf-8")
    return pd.Index(genes)


def load_visium_sample(sample_dir: Path) -> ad.AnnData:
    outs = sample_dir / "outs"
    adata = sc.read_visium(outs)
    adata.var_names_make_unique()
    adata = filter_mt_genes(adata)
    meta = parse_sample_metadata(sample_dir.name)
    for k, v in meta.items():
        adata.obs[k] = v
    old_key = next(iter(adata.uns["spatial"].keys()))
    adata.uns["spatial"] = {meta["sample_id"]: adata.uns["spatial"][old_key]}
    return adata


def run_one_sample(
    sample_dir: Path,
    ref_sig: pd.DataFrame,
    selected_genes: pd.Index,
    output_dir: Path,
    accelerator: str,
    max_epochs: int,
    batch_size: int,
) -> dict[str, object]:
    sample_id = sample_dir.name.lstrip("#")
    sample_out = output_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)
    out_h5ad = sample_out / "cell2loc.h5ad"
    out_csv = sample_out / "cell2loc_q05.csv"

    if out_h5ad.exists() and out_csv.exists():
        q05 = pd.read_csv(out_csv, index_col=0)
        return {
            "sample_id": sample_id,
            "n_spots": int(q05.shape[0]),
            "n_cell_types": int(q05.shape[1]),
            "status": "skipped_existing",
        }

    print(f"[sample] {sample_id}")
    adata = load_visium_sample(sample_dir)
    shared = adata.var_names.intersection(ref_sig.index).intersection(selected_genes)
    adata = adata[:, shared].copy()
    sig = ref_sig.loc[shared, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata, batch_key=None)
    model = cell2location.models.Cell2location(
        adata,
        cell_state_df=sig,
        N_cells_per_location=8,
        detection_alpha=20,
    )
    model.train(
        max_epochs=max_epochs,
        accelerator=accelerator,
        batch_size=batch_size,
        train_size=1,
    )
    model.save(sample_out / "model", overwrite=True)

    adata = model.export_posterior(
        adata,
        sample_kwargs={"num_samples": 1000, "batch_size": min(batch_size, adata.n_obs)},
    )
    adata.write(out_h5ad, compression="gzip")

    cell_types = list(adata.uns["mod"]["factor_names"])
    q05 = pd.DataFrame(adata.obsm["q05_cell_abundance_w_sf"], index=adata.obs_names, columns=cell_types)
    means = pd.DataFrame(adata.obsm["means_cell_abundance_w_sf"], index=adata.obs_names, columns=cell_types)
    q05.to_csv(out_csv)
    means.to_csv(sample_out / "cell2loc_mean.csv")
    adata.obs.to_csv(sample_out / "obs.csv")

    with open(sample_out / "meta.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "sample_id": sample_id,
                "n_spots": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "cell_types": cell_types,
                "max_epochs": max_epochs,
                "batch_size": batch_size,
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )
    return {
        "sample_id": sample_id,
        "n_spots": int(adata.n_obs),
        "n_cell_types": len(cell_types),
        "status": "done",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run cell2location on Visium 2 samples one by one.")
    parser.add_argument("--visium-root", required=True, type=Path)
    parser.add_argument("--ref-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample-prefix", default="#UKF")
    parser.add_argument("--max-epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    torch.set_float32_matmul_precision("medium")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    args.output.mkdir(parents=True, exist_ok=True)
    ref_sig = load_reference_signatures(args.ref_dir)
    selected_genes = load_selected_genes(args.ref_dir)

    sample_dirs = sorted(
        p for p in args.visium_root.iterdir() if p.is_dir() and p.name.startswith(args.sample_prefix)
    )
    if args.limit > 0:
        sample_dirs = sample_dirs[: args.limit]

    rows = []
    for sample_dir in sample_dirs:
        rows.append(
            run_one_sample(
                sample_dir=sample_dir,
                ref_sig=ref_sig,
                selected_genes=selected_genes,
                output_dir=args.output,
                accelerator=accelerator,
                max_epochs=args.max_epochs,
                batch_size=args.batch_size,
            )
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output / "summary.csv", index=False)
    with open(args.output / "run.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "visium_root": str(args.visium_root),
                "ref_dir": str(args.ref_dir),
                "accelerator": accelerator,
                "selected_genes": int(len(selected_genes)),
                "max_epochs": args.max_epochs,
                "batch_size": args.batch_size,
                "n_samples": len(rows),
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )
    print("[done]", len(rows))


if __name__ == "__main__":
    main()
