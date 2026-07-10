from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_bbknn_analysis" / "gbm_merged_raw_counts_ilc_qc_bbknn_clustered.h5ad"
OUTDIR = ROOT / "gbm_bbknn_analysis" / "featureplots_major_groups"

MARKER_SETS = {
    "Glial": ["MAP1B", "NCAM1", "GPM6B"],
    "Mac": ["SPP1", "HLA-DRA", "P2RY12", "CD163"],
    "Monocytes": ["FCN1", "VCAN", "LYZ"],
    "Neutrophils": ["S100A8", "S100A9", "FCGR3B", "DEFA1"],
    "T_cells": ["CD3D", "TRAC", "IL7R", "CCL5"],
    "NK_ILC": ["NKG7", "GNLY", "KLRD1", "KLRB1"],
    "B_cells": ["MS4A1", "CD79A", "CD79B"],
    "ILC1": ["IL7R", "TBX21", "CXCR6", "XCL1", "XCL2", "IFNG"],
    "ILC2": ["IL7R", "KLRB1", "GATA3", "RORA", "IL1RL1", "PTGDR2", "AREG", "ICOS"],
    "ILC3": ["IL7R", "KLRB1", "RORC", "AHR", "IL23R", "KIT", "NCR2", "LTB"],
}


def save_plot(adata, genes: list[str], output: Path) -> None:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return
    sc.pl.umap(
        adata,
        color=present,
        cmap="RdBu_r",
        vmin=0,
        vmax="p99",
        show=False,
    )
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    for group, genes in MARKER_SETS.items():
        save_plot(adata, genes, OUTDIR / f"{group}_featureplot.png")

    all_genes = []
    for genes in MARKER_SETS.values():
        all_genes.extend(genes[:2])
    deduped = []
    seen = set()
    for gene in all_genes:
        if gene not in seen:
            deduped.append(gene)
            seen.add(gene)
    save_plot(adata, deduped, OUTDIR / "all_major_groups_featureplot.png")
    print(f"written {OUTDIR}")


if __name__ == "__main__":
    main()
