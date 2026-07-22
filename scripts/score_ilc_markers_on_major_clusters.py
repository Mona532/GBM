from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered.h5ad"
OUTDIR = ROOT / "gbm_bbknn_analysis" / "ilc_scoring"
OUTPUT = OUTDIR / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered_ilc_scored.h5ad"

MARKERS = {
    "pan_ilc": ["IL7R", "KLRB1", "KIT", "RORA", "TCF7", "AHR", "ICOS", "CXCR6"],
    "ilc1_like": ["IL7R", "TBX21", "CXCR6", "XCL1", "XCL2", "IFNG", "TNFRSF18", "TYROBP"],
    "ilc2_like": ["IL7R", "KLRB1", "GATA3", "RORA", "IL1RL1", "PTGDR2", "AREG", "ICOS"],
    "ilc3_like": ["IL7R", "KLRB1", "RORC", "AHR", "IL23R", "KIT", "NCR2", "LTB"],
}


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


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT)
    if "lognorm" not in adata.layers:
        raise RuntimeError("Expected lognorm layer in clustered H5AD")

    work = adata.copy()
    work.X = work.layers["lognorm"].copy()

    marker_rows = []
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in work.var_names]
        missing = [g for g in genes if g not in work.var_names]
        if len(present) < 2:
            raise RuntimeError(f"Too few markers present for {name}: {present}")
        sc.tl.score_genes(work, gene_list=present, score_name=f"{name}_score", use_raw=False)
        marker_rows.append(
            {
                "signature": name,
                "n_markers_total": len(genes),
                "n_markers_present": len(present),
                "markers_present": ",".join(present),
                "markers_missing": ",".join(missing),
            }
        )

    marker_df = pd.DataFrame(marker_rows)
    marker_df.to_csv(OUTDIR / "ilc_marker_sets.tsv", sep="\t", index=False)

    score_cols = [f"{name}_score" for name in MARKERS]
    cluster_summary = (
        work.obs.groupby("seurat_clusters", observed=True)[score_cols]
        .agg(["mean", "median"])
        .sort_index()
    )
    cluster_summary.to_csv(OUTDIR / "ilc_cluster_score_summary.tsv", sep="\t")

    mean_scores = work.obs.groupby("seurat_clusters", observed=True)[score_cols].mean().sort_index()
    mean_scores.to_csv(OUTDIR / "ilc_cluster_mean_scores.tsv", sep="\t")

    zscores = mean_scores.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) if col.std(ddof=0) != 0 else 1), axis=0)
    zscores.to_csv(OUTDIR / "ilc_cluster_mean_scores_zscore.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(mean_scores, cmap="viridis", ax=ax)
    ax.set_title("ILC scores")
    fig.tight_layout()
    fig.savefig(OUTDIR / "ilc_cluster_mean_scores_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(zscores, cmap="coolwarm", center=0, ax=ax)
    ax.set_title("ILC scores z")
    fig.tight_layout()
    fig.savefig(OUTDIR / "ilc_cluster_mean_scores_zscore_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    winner_rows = []
    for sig in MARKERS:
        col = f"{sig}_score"
        best_cluster = mean_scores[col].idxmax()
        winner_rows.append(
            {
                "signature": sig,
                "top_cluster": str(best_cluster),
                "top_mean_score": float(mean_scores.loc[best_cluster, col]),
                "top_zscore": float(zscores.loc[best_cluster, col]),
            }
        )
    pd.DataFrame(winner_rows).to_csv(OUTDIR / "ilc_top_clusters.tsv", sep="\t", index=False)

    for col in score_cols:
        adata.obs[col] = work.obs[col].values
    adata = sanitize_for_h5ad_write(adata)
    adata.write_h5ad(OUTPUT, compression="gzip")
    print(f"written {OUTPUT}")
    print("top clusters")
    print(pd.DataFrame(winner_rows).to_string(index=False))


if __name__ == "__main__":
    main()
