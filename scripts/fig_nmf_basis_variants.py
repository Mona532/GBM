from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    }
)

ROOT = Path(r"E:/GBM/results")
ECO_ORDER = ["E1", "E2", "E3", "E4"]
ECO_SHORT = {
    "E1": "Lymphocyte",
    "E2": "ILC-enriched",
    "E3": "Myeloid-vascular",
    "E4": "Glial-CD4",
}
ECO_COLORS = {
    "E1": "#c44e52",
    "E2": "#4c72b0",
    "E3": "#55a868",
    "E4": "#dd8452",
}
CELLTYPE_ORDER = [
    "B",
    "Plasma",
    "CD4_T",
    "CD8_T",
    "NK",
    "Dendritic",
    "ILC1",
    "ILC2",
    "ILC3",
    "Macrophage",
    "Vascular",
    "Glial",
    "Glioma",
]


def save_multi(fig, stem: str):
    fig.savefig(ROOT / f"{stem}.jpg", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def draw_heatmap(mat, row_labels, col_labels, title, stem, cmap, vmin=None, vmax=None, cbar_label="Value"):
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=25, ha="right", fontsize=8)
    for idx, eco in enumerate(ECO_ORDER):
        ax.get_xticklabels()[idx].set_color(ECO_COLORS[eco])
    ax.set_title(title, fontsize=10, pad=8)
    cbar = plt.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
    cbar.set_label(cbar_label, fontsize=8)
    fig.tight_layout()
    save_multi(fig, stem)


def main():
    basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
    basis = basis.loc[CELLTYPE_ORDER, ECO_ORDER]

    row_scaled = StandardScaler().fit_transform(basis.values.T).T
    row_v = np.abs(row_scaled).max()
    draw_heatmap(
        row_scaled,
        CELLTYPE_ORDER,
        [ECO_SHORT[e] for e in ECO_ORDER],
        "NMF Basis: Row-scaled",
        "fig_nmf_basis_row_scaled",
        "RdBu_r",
        vmin=-row_v,
        vmax=row_v,
        cbar_label="Row z-score",
    )

    col_scaled = StandardScaler().fit_transform(basis.values)
    col_v = np.abs(col_scaled).max()
    draw_heatmap(
        col_scaled,
        CELLTYPE_ORDER,
        [ECO_SHORT[e] for e in ECO_ORDER],
        "NMF Basis: Column-scaled",
        "fig_nmf_basis_col_scaled",
        "RdBu_r",
        vmin=-col_v,
        vmax=col_v,
        cbar_label="Column z-score",
    )

    colsum_norm = basis.div(basis.sum(axis=0), axis=1)
    draw_heatmap(
        colsum_norm.values,
        CELLTYPE_ORDER,
        [ECO_SHORT[e] for e in ECO_ORDER],
        "NMF Basis: Column-normalized Composition",
        "fig_nmf_basis_colsum_norm",
        "YlOrRd",
        vmin=0,
        vmax=float(np.quantile(colsum_norm.values, 0.98)),
        cbar_label="Composition weight",
    )

    global_log = np.log1p(basis.values)
    draw_heatmap(
        global_log,
        CELLTYPE_ORDER,
        [ECO_SHORT[e] for e in ECO_ORDER],
        "NMF Basis: Global log1p",
        "fig_nmf_basis_global_log1p",
        "YlOrRd",
        vmin=float(global_log.min()),
        vmax=float(np.quantile(global_log, 0.98)),
        cbar_label="log1p(weight)",
    )

    print("Saved NMF basis variants.")


if __name__ == "__main__":
    main()
