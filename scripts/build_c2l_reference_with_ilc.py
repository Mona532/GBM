from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import concat as ad_concat
from cell2location.utils.filtering import filter_genes


ILC_MARKERS = [
    "TBX21",
    "XCL1",
    "XCL2",
    "IFNG",
    "NCR1",
    "GATA3",
    "IL7R",
    "RORA",
    "AREG",
    "IL1RL1",
    "KLRB1",
    "RORC",
    "KIT",
    "IL23R",
    "CCR6",
    "AHR",
    "LTB",
]

LYMPHOID_FINE_MAP = {
    "NK cells 1": "NK",
    "NK cells 2": "NK",
    "CD8+ T cells": "CD8_T",
    "CD8+ T cells (cytotoxic)": "CD8_T",
    "HSP-response T cells": "CD8_T",
    "IFN-response T cells": "CD8_T",
    "Naïve T cells": "CD8_T",
    "Nave T cells": "CD8_T",
    "B cells": "B",
    "Dendritic cells": "Dendritic",
}

KEEP_COARSE = {
    "AC-gliosis-like",
    "AC-progenitor-like",
    "Astrocytes",
    "Gliosis-like",
    "Hypoxic",
    "Myeloid",
    "NPC-neuronal-like",
    "Neurons (Exc)",
    "Neurons (Inh)",
    "OPC-NPC-like",
    "OPC-like",
    "OPC-neuronal-like",
    "OPCs",
    "Oligodendrocytes",
    "Proliferative",
    "Vascular-associated",
}


def ensure_raw_counts(adata: ad.AnnData, name: str) -> ad.AnnData:
    adata = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    x = adata.X
    vals = x.data[:10000] if sp.issparse(x) else np.ravel(x)[:10000]
    vals = np.asarray(vals)
    if vals.size and np.mean(np.abs(vals - np.round(vals)) > 1e-6) > 0:
        raise ValueError(f"{name} does not look like raw count data")
    return adata


def set_main_ref_gene_symbols(adata: ad.AnnData) -> ad.AnnData:
    if "SYMBOL" in adata.var.columns:
        symbols = adata.var["SYMBOL"].astype(str).to_numpy()
        fallback = adata.var_names.astype(str).to_numpy()
        names = np.where((symbols != "nan") & (symbols != ""), symbols, fallback)
        adata.var_names = pd.Index(names)
    adata.var_names_make_unique()
    return adata


def filter_mt_genes(adata: ad.AnnData) -> ad.AnnData:
    keep = ~adata.var_names.str.upper().str.startswith("MT-")
    return adata[:, keep].copy()


def choose_cells_by_label(
    labels: pd.Series,
    keep_labels: set[str],
    max_cells_per_label: int,
    seed: int,
) -> np.ndarray:
    labels = labels.astype(str)
    keep_mask = labels.isin(keep_labels)
    kept = labels.loc[keep_mask]
    rng = np.random.default_rng(seed)
    keep_idx: list[np.ndarray] = []
    for _, idx in kept.groupby(kept).groups.items():
        idx = np.asarray(list(idx))
        if max_cells_per_label > 0 and len(idx) > max_cells_per_label:
            idx = rng.choice(idx, size=max_cells_per_label, replace=False)
        keep_idx.append(idx)
    return np.concatenate(keep_idx)


def build_main_labels(adata: ad.AnnData) -> tuple[ad.AnnData, set[str]]:
    coarse = adata.obs["annotation_coarse"].astype(str)
    fine = adata.obs["annotation_granular"].astype(str)
    labels = coarse.copy()
    fine_hits = fine.map(LYMPHOID_FINE_MAP)
    labels.loc[fine_hits.notna()] = fine_hits.loc[fine_hits.notna()]
    keep_labels = KEEP_COARSE | set(LYMPHOID_FINE_MAP.values())
    keep_idx = choose_cells_by_label(labels, keep_labels, 0, 1234)
    adata = adata[keep_idx].copy()
    adata.obs["c2l_label"] = labels.loc[keep_idx].to_numpy()
    adata.obs["c2l_batch"] = adata.obs["donor_id"].astype(str)
    adata.obs["ref_source"] = "GBM_space_snRNA"
    return adata, keep_labels


def build_ilc_labels(adata: ad.AnnData) -> ad.AnnData:
    labels = adata.obs["ilc_subtype"].astype(str)
    keep = labels.isin({"ILC1", "ILC2", "ILC3", "NK"})
    adata = adata[keep.to_numpy()].copy()
    adata.obs["c2l_label"] = labels.loc[keep].to_numpy()
    adata.obs["c2l_batch"] = "ILC_" + adata.obs["sample"].astype(str)
    adata.obs["ref_source"] = "GBM_ilc"
    return adata


def write_summary(adata: ad.AnnData, out_dir: Path) -> None:
    summary = (
        adata.obs.groupby(["c2l_label", "ref_source"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["c2l_label", "ref_source"])
    )
    summary.to_csv(out_dir / "reference_label_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a full cell2location reference with ILC-enhanced lymphoid labels.")
    parser.add_argument("--main-ref", required=True, type=Path)
    parser.add_argument("--ilc-ref", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-cells-per-label", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--min-cells-per-gene", type=int, default=3)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    main_ref = ad.read_h5ad(args.main_ref)
    main_ref = ensure_raw_counts(main_ref, "main_ref")
    main_ref = set_main_ref_gene_symbols(main_ref)
    main_ref = filter_mt_genes(main_ref)
    coarse = main_ref.obs["annotation_coarse"].astype(str)
    fine = main_ref.obs["annotation_granular"].astype(str)
    main_labels = coarse.copy()
    fine_hits = fine.map(LYMPHOID_FINE_MAP)
    main_labels.loc[fine_hits.notna()] = fine_hits.loc[fine_hits.notna()]
    main_keep_labels = KEEP_COARSE | set(LYMPHOID_FINE_MAP.values())
    main_keep_idx = choose_cells_by_label(
        main_labels,
        main_keep_labels,
        args.max_cells_per_label,
        args.seed,
    )
    main_ref = main_ref[main_keep_idx].copy()
    main_ref.obs["c2l_label"] = main_labels.loc[main_keep_idx].to_numpy()
    main_ref.obs["c2l_batch"] = main_ref.obs["donor_id"].astype(str)
    main_ref.obs["ref_source"] = "GBM_space_snRNA"

    ilc_ref = ad.read_h5ad(args.ilc_ref)
    ilc_ref = ensure_raw_counts(ilc_ref, "ilc_ref")
    ilc_ref.var_names_make_unique()
    ilc_ref = filter_mt_genes(ilc_ref)
    ilc_ref = build_ilc_labels(ilc_ref)

    shared = main_ref.var_names.intersection(ilc_ref.var_names)
    if len(shared) < 2000:
        raise RuntimeError(f"Too few shared genes between references: {len(shared)}")

    main_ref = main_ref[:, shared].copy()
    ilc_ref = ilc_ref[:, shared].copy()

    combined = ad_concat(
        {"main": main_ref, "ilc": ilc_ref},
        axis=0,
        join="inner",
        merge="same",
        label="concat_source",
        index_unique="-",
    )
    if "concat_source" in combined.obs:
        del combined.obs["concat_source"]

    sc.pp.filter_genes(combined, min_cells=args.min_cells_per_gene)

    selected = pd.Index(
        filter_genes(
            combined,
            cell_count_cutoff=15,
            cell_percentage_cutoff2=0.05,
            nonz_mean_cutoff=1.12,
        )
    )
    forced = pd.Index([g for g in ILC_MARKERS if g in combined.var_names])
    selected = selected.union(forced)

    combined.write(args.output / "reference_with_ilc.h5ad", compression="gzip")
    (args.output / "selected_genes_ilc.txt").write_text("\n".join(map(str, selected)), encoding="utf-8")
    write_summary(combined, args.output)

    meta = {
        "main_ref": str(args.main_ref),
        "ilc_ref": str(args.ilc_ref),
        "output_h5ad": str(args.output / "reference_with_ilc.h5ad"),
        "selected_genes": int(len(selected)),
        "forced_ilc_markers_kept": [g for g in forced],
        "max_cells_per_label": args.max_cells_per_label,
        "seed": args.seed,
        "n_cells": int(combined.n_obs),
        "n_genes": int(combined.n_vars),
        "labels": combined.obs["c2l_label"].astype(str).value_counts().to_dict(),
    }
    with open(args.output / "build_meta.json", "w", encoding="utf-8") as handle:
        json.dump(meta, handle, ensure_ascii=False, indent=2)

    print("[done]")


if __name__ == "__main__":
    main()
