from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_lymphoid_fine_analysis" / "gbm_lymphoid_fine_bbknn_clustered.h5ad"
OUTDIR = ROOT / "gbm_lymphoid_fine_analysis"
PLOT_OUT = OUTDIR / "mt_percent_by_lymphoid_fine_clusters.png"
SUMMARY_OUT = OUTDIR / "mt_percent_by_lymphoid_fine_clusters_summary.tsv"


def main() -> None:
    adata = sc.read_h5ad(INPUT)

    if "percent.mt" in adata.obs.columns:
        mt_col = "percent.mt"
    elif "pct_counts_mt" in adata.obs.columns:
        mt_col = "pct_counts_mt"
    else:
        raise RuntimeError("No mitochondrial percentage column found in obs")

    cluster_col = "lymphoid_fine_clusters"
    if cluster_col not in adata.obs.columns:
        raise RuntimeError(f"Missing {cluster_col} in obs")

    plot_df = adata.obs[[cluster_col, mt_col]].copy()
    plot_df[cluster_col] = plot_df[cluster_col].astype(str)
    plot_df[mt_col] = pd.to_numeric(plot_df[mt_col], errors="coerce")

    summary = (
        plot_df.groupby(cluster_col, observed=True)[mt_col]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
        .sort_values(cluster_col)
    )
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.violinplot(data=plot_df, x=cluster_col, y=mt_col, inner="box", cut=0, ax=ax, color="#72b7b2")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("MT%")
    ax.set_title("MT%")
    fig.tight_layout()
    fig.savefig(PLOT_OUT, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"written {PLOT_OUT}")
    print(f"written {SUMMARY_OUT}")


if __name__ == "__main__":
    main()
