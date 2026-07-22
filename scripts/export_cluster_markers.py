from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered.h5ad"
OUTDIR = ROOT / "gbm_bbknn_analysis"
ALL_OUT = OUTDIR / "cluster_markers_all.tsv"
TOP_OUT = OUTDIR / "cluster_markers_top20.tsv"


def main() -> None:
    adata = sc.read_h5ad(INPUT)
    if "rank_genes_groups" not in adata.uns:
        raise RuntimeError("Missing rank_genes_groups in clustered H5AD")

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
    all_markers.to_csv(ALL_OUT, sep="\t", index=False)

    top_markers = (
        all_markers.sort_values(["cluster", "rank"], ascending=[True, True])
        .groupby("cluster", as_index=False)
        .head(20)
        .reset_index(drop=True)
    )
    top_markers.to_csv(TOP_OUT, sep="\t", index=False)

    print(f"written {ALL_OUT}")
    print(f"written {TOP_OUT}")
    print(f"clusters {len(groups)}")


if __name__ == "__main__":
    main()
