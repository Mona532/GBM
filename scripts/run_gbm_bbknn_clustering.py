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
INPUT = ROOT / "gbm_merged_raw_counts_ilc_qc_symbol_safe.h5ad"
OUTDIR = ROOT / "gbm_bbknn_analysis"
OUTPUT = OUTDIR / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered.h5ad"


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


def save_violin(adata: ad.AnnData, keys: list[str], output: Path, groupby: str | None = None) -> None:
    sc.pl.violin(adata, keys=keys, groupby=groupby, show=False, multi_panel=True)
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_umap(adata: ad.AnnData, color, output: Path, title: str | None = None, wspace: float | None = None) -> None:
    sc.pl.umap(adata, color=color, show=False, legend_fontsize=8, frameon=False, wspace=wspace)
    fig = plt.gcf()
    if title:
        fig.suptitle(title)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_pca_variance(adata: ad.AnnData, output: Path) -> None:
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_hvg_plot(adata: ad.AnnData, output: Path) -> None:
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_source_barplot(adata: ad.AnnData, output: Path) -> None:
    counts = adata.obs["source_for_bbknn"].astype(str).value_counts()
    fig, ax = plt.subplots(figsize=(7, 4))
    counts.plot(kind="bar", ax=ax, color="#4c78a8")
    ax.set_ylabel("Cells")
    ax.set_xlabel("Source")
    ax.set_title("Source")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_qc_post_plot(adata: ad.AnnData, output: Path) -> None:
    cols = ["nCount_RNA", "nFeature_RNA"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, col in zip(axes, cols):
        ax.hist(pd.to_numeric(adata.obs[col], errors="coerce").fillna(0), bins=50, color="#72b7b2")
        ax.set_title(col)
        ax.set_xlabel(col)
        ax.set_ylabel("Cells")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_qc_retention_plot(adata: ad.AnnData, output: Path) -> None:
    current = adata.obs["dataset"].astype(str).value_counts().reindex(["GBM-ILC1", "GBM-ILC2", "GBM01", "GBM02"]).fillna(0)
    original = pd.Series({"GBM-ILC1": 3000, "GBM-ILC2": 2284, "GBM01": 6199, "GBM02": 9759})
    table = pd.DataFrame({"before_qc": original, "after_qc_or_unchanged": current.astype(int)})
    fig, ax = plt.subplots(figsize=(8, 4))
    table.plot(kind="bar", ax=ax, color=["#bab0ab", "#e15759"])
    ax.set_ylabel("Cells")
    ax.set_title("QC")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    source_col = "Source" if "Source" in adata.obs.columns else ("source" if "source" in adata.obs.columns else None)
    sample_str = adata.obs["sample"].astype(str)
    if source_col is None:
        adata.obs["source_for_bbknn"] = sample_str
    else:
        src = adata.obs[source_col].astype(str).replace({"": np.nan, "nan": np.nan, "None": np.nan})
        adata.obs["source_for_bbknn"] = np.where(
            src.notna(),
            sample_str + "_" + src.astype(str),
            sample_str,
        )
        adata.obs["source_for_bbknn"] = adata.obs["source_for_bbknn"].astype(str)

    save_qc_retention_plot(adata, OUTDIR / "step0_qc_retention.png")
    save_qc_post_plot(adata, OUTDIR / "step0_qc_post_metrics.png")
    save_source_barplot(adata, OUTDIR / "step1_source_for_bbknn_counts.png")

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    save_violin(adata, ["nCount_RNA", "nFeature_RNA"], OUTDIR / "step2_lognorm_violin_by_dataset.png", groupby="dataset")

    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="source_for_bbknn")
    save_hvg_plot(adata, OUTDIR / "step3_hvg_3000.png")

    adata.uns["log1p"] = {"base": None}
    sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    save_pca_variance(adata, OUTDIR / "step4_pca_variance_ratio.png")

    bbknn.bbknn(adata, batch_key="source_for_bbknn", use_rep="X_pca")
    sc.tl.umap(adata)
    save_umap(adata, "source_for_bbknn", OUTDIR / "step5_umap_source_for_bbknn.png", "Source")
    save_umap(adata, "dataset", OUTDIR / "step5_umap_dataset.png", "Dataset")

    sc.tl.leiden(adata, resolution=1.0, key_added="seurat_clusters")
    save_umap(adata, ["seurat_clusters"], OUTDIR / "step6_umap_leiden_res1.0.png", "Leiden")

    sc.tl.rank_genes_groups(adata, groupby="seurat_clusters", method="wilcoxon", n_genes=50)
    sc.pl.rank_genes_groups(adata, n_genes=20, show=False)
    fig = plt.gcf()
    fig.savefig(OUTDIR / "step7_rank_genes_groups_top20.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    adata = sanitize_for_h5ad_write(adata)
    adata.write_h5ad(OUTPUT, compression="gzip")
    print(f"written {OUTPUT}")
    print(f"shape {adata.shape}")


if __name__ == "__main__":
    main()
