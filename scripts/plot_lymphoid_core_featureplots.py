from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "gbm_lymphoid_fine_analysis" / "gbm_lymphoid_fine_bbknn_clustered.h5ad"
OUTDIR = ROOT / "gbm_lymphoid_fine_analysis" / "featureplots_core"

MARKER_SETS = {
    "ILC_core": [
        "TBX21", "IFNG", "TNF", "CXCR3", "IL12RB2", "IL18R1", "XCL1", "XCL2",
        "NCR1", "NCR3", "RUNX3", "HOPX", "GATA3", "RORA", "PTGDR2", "HPGDS",
        "IL1RL1", "IL17RB", "AREG", "RORC", "AHR", "IL23R", "KIT", "CCR6", "IL22",
    ],
    "NK_core": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "KLRD1", "KLRF1", "FCGR3A", "EOMES"],
    "T_core": ["CD3D", "CD3E", "TRAC", "IL7R", "LTB", "ITK", "CCL5", "BACH2", "RORA"],
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
        frameon=False,
    )
    fig = plt.gcf()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    for name, genes in MARKER_SETS.items():
        save_plot(adata, genes, OUTDIR / f"{name}_featureplot.png")

    print(f"written {OUTDIR}")


if __name__ == "__main__":
    main()
