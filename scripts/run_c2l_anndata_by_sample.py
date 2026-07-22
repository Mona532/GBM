from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import cell2location
import numpy as np
import pandas as pd
import torch


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

    from cell2location.utils.filtering import filter_genes

    ref_posterior = ad.read_h5ad(ref_dir / "reference_posterior.h5ad")
    genes = filter_genes(
        ref_posterior,
        cell_count_cutoff=15,
        cell_percentage_cutoff2=0.05,
        nonz_mean_cutoff=1.12,
    )
    genes_txt.write_text("\n".join(map(str, genes)), encoding="utf-8")
    return pd.Index(genes)


def load_sample_h5ad(h5ad_path: Path) -> ad.AnnData:
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names_make_unique()
    mask = adata.var["feature_types"] == "Gene Expression"
    adata = adata[:, mask].copy()
    return adata


def run_one_sample(
    h5ad_path: Path,
    ref_sig: pd.DataFrame,
    selected_genes: pd.Index,
    output_dir: Path,
    accelerator: str,
    max_epochs: int,
    batch_size: int,
) -> dict[str, object]:
    sample_id = h5ad_path.stem
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
    adata = load_sample_h5ad(h5ad_path)

    shared = adata.var_names.intersection(ref_sig.index).intersection(selected_genes)
    print(f"  genes shared with reference: {len(shared)}")
    adata = adata[:, shared].copy()
    sig = ref_sig.loc[shared, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata, batch_key=None)
    model = cell2location.models.Cell2location(
        adata,
        cell_state_df=sig,
        N_cells_per_location=8,
        detection_alpha=200,
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
    q05_arr = adata.obsm["q05_cell_abundance_w_sf"]
    means_arr = adata.obsm["means_cell_abundance_w_sf"]
    if hasattr(q05_arr, "values"):
        q05_arr = q05_arr.values
    if hasattr(means_arr, "values"):
        means_arr = means_arr.values
    q05 = pd.DataFrame(q05_arr, index=adata.obs_names, columns=cell_types)
    means = pd.DataFrame(means_arr, index=adata.obs_names, columns=cell_types)
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
    print(f"  done: {adata.n_obs} spots, {len(cell_types)} cell types")
    return {
        "sample_id": sample_id,
        "n_spots": int(adata.n_obs),
        "n_cell_types": len(cell_types),
        "status": "done",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run cell2location on pre-processed .h5ad Visium samples one by one."
    )
    parser.add_argument("--anndata-dir", required=True, type=Path)
    parser.add_argument("--ref-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-epochs", type=int, default=2000)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    torch.set_float32_matmul_precision("medium")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    args.output.mkdir(parents=True, exist_ok=True)

    ref_sig = load_reference_signatures(args.ref_dir)
    selected_genes = load_selected_genes(args.ref_dir)

    h5ad_paths = sorted(args.anndata_dir.glob("*.h5ad"))
    if args.limit > 0:
        h5ad_paths = h5ad_paths[: args.limit]

    print(f"Samples to process: {len(h5ad_paths)}")
    print(f"Reference cell types: {ref_sig.shape[1]}")
    print(f"Selected genes: {len(selected_genes)}")
    print(f"Accelerator: {accelerator}")

    rows = []
    for h5ad_path in h5ad_paths:
        rows.append(
            run_one_sample(
                h5ad_path=h5ad_path,
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
                "anndata_dir": str(args.anndata_dir),
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
    done = sum(1 for r in rows if r["status"] == "done")
    skipped = sum(1 for r in rows if r["status"] == "skipped_existing")
    print(f"[done] {done} processed, {skipped} skipped (already existed)")


if __name__ == "__main__":
    main()
