from __future__ import annotations

from pathlib import Path

import anndata as ad
import bbknn
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered_symbolsafe_lognorm_for_celltypist.h5ad"
OUTDIR = ROOT / "gbm_lymphoid_fine_analysis"
OUTPUT = OUTDIR / "gbm_lymphoid_fine_bbknn_clustered.h5ad"
ALL_MARKERS = OUTDIR / "cluster_markers_all.tsv"
TOP_MARKERS = OUTDIR / "cluster_markers_top30.tsv"
CLUSTERS_KEEP = {"3", "8", "9"}


def sanitize_for_h5ad_write(adata: ad.AnnData) -> ad.AnnData:
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = True
    adata.obs_names = pd.Index(adata.obs_names.astype(str))
    adata.var_names = pd.Index(adata.var_names.astype(str))
    for frame in [adata.obs, adata.var]:
        frame.index = pd.Index(frame.index.astype(str))
        for col in frame.columns:
            if isinstance(frame[col].dtype, pd.CategoricalDtype):
                continue
            if pd.api.types.is_string_dtype(frame[col]):
                frame[col] = frame[col].fillna("").astype(str).astype(object)
    return adata


def save_umap(adata: ad.AnnData, color, output: Path, title: str) -> None:
    sc.pl.umap(adata, color=color, show=False, frameon=False, legend_fontsize=8)
    fig = plt.gcf()
    fig.suptitle(title)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_pca_variance(adata: ad.AnnData, output: Path) -> None:
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=30, show=False)
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_hvg_plot(adata: ad.AnnData, output: Path) -> None:
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def export_markers(adata: ad.AnnData) -> None:
    rg = adata.uns["rank_genes_groups"]
    groups = list(rg["names"].dtype.names)
    rows = []
    for group in groups:
        names = rg["names"][group]
        scores = rg["scores"][group] if "scores" in rg else [None] * len(names)
        pvals = rg["pvals"][group] if "pvals" in rg else [None] * len(names)
        pvals_adj = rg["pvals_adj"][group] if "pvals_adj" in rg else [None] * len(names)
        logfcs = rg["logfoldchanges"][group] if "logfoldchanges" in rg else [None] * len(names)
        for rank, gene in enumerate(names, start=1):
            rows.append(
                {
                    "cluster": str(group),
                    "rank": rank,
                    "gene": gene,
                    "score": float(scores[rank - 1]) if pd.notna(scores[rank - 1]) else None,
                    "logfoldchange": float(logfcs[rank - 1]) if pd.notna(logfcs[rank - 1]) else None,
                    "pval": float(pvals[rank - 1]) if pd.notna(pvals[rank - 1]) else None,
                    "pval_adj": float(pvals_adj[rank - 1]) if pd.notna(pvals_adj[rank - 1]) else None,
                }
            )
    all_markers = pd.DataFrame(rows)
    symbol_markers = all_markers.loc[~all_markers["gene"].astype(str).str.startswith("ENSG")].copy()
    symbol_markers.to_csv(ALL_MARKERS, sep="\t", index=False)
    (
        symbol_markers.sort_values(["cluster", "rank"], ascending=[True, True])
        .groupby("cluster", as_index=False)
        .head(30)
        .reset_index(drop=True)
        .to_csv(TOP_MARKERS, sep="\t", index=False)
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT)
    adata.var["gene_symbol"] = adata.var_names.astype(str)
    keep_genes = ~adata.var["gene_symbol"].astype(str).str.startswith("ENSG")
    adata = adata[:, keep_genes].copy()
    keep = adata.obs["seurat_clusters"].astype(str).isin(CLUSTERS_KEEP)
    adata = adata[keep].copy()

    adata.obs["lymphoid_parent"] = np.where(
        adata.obs["seurat_clusters"].astype(str).eq("8"),
        "NK_ILC",
        "T_cells",
    )

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sample_str = adata.obs["sample"].astype(str)
    source_col = "Source" if "Source" in adata.obs.columns else ("source" if "source" in adata.obs.columns else None)
    if source_col is None:
        adata.obs["source_for_bbknn"] = sample_str
    else:
        src = adata.obs[source_col].astype(str).replace({"": np.nan, "nan": np.nan, "None": np.nan})
        adata.obs["source_for_bbknn"] = np.where(src.notna(), sample_str + "_" + src.astype(str), sample_str).astype(str)

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    save_hvg_plot(adata, OUTDIR / "step1_hvg_3000.png")

    adata.uns["log1p"] = {"base": None}
    sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    save_pca_variance(adata, OUTDIR / "step2_pca_variance_ratio.png")

    bbknn.bbknn(adata, batch_key="source_for_bbknn", use_rep="X_pca")
    sc.tl.umap(adata)
    save_umap(adata, ["lymphoid_parent", "dataset"], OUTDIR / "step3_umap_parent_dataset.png", "UMAP")

    sc.tl.leiden(adata, resolution=1.5, key_added="lymphoid_fine_clusters")
    save_umap(adata, ["lymphoid_fine_clusters"], OUTDIR / "step4_umap_leiden_res1.5.png", "Leiden")

    sc.tl.rank_genes_groups(adata, groupby="lymphoid_fine_clusters", method="wilcoxon", n_genes=80)
    sc.pl.rank_genes_groups(adata, n_genes=25, show=False)
    fig = plt.gcf()
    fig.savefig(OUTDIR / "step5_rank_genes_groups_top25.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    export_markers(adata)
    adata = sanitize_for_h5ad_write(adata)
    adata.write_h5ad(OUTPUT, compression="gzip")
    print(f"written {OUTPUT}")
    print(f"shape {adata.shape}")


if __name__ == "__main__":
    main()
