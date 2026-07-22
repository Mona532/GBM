from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
        "legend.frameon": False,
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
CELLTYPE_COLORS = {
    "B": "#4e79a7",
    "Plasma": "#a0cbe8",
    "CD4_T": "#e15759",
    "CD8_T": "#f28e2b",
    "NK": "#ffbe7d",
    "Dendritic": "#76b7b2",
    "ILC1": "#59a14f",
    "ILC2": "#8cd17d",
    "ILC3": "#b6992d",
    "Macrophage": "#9c755f",
    "Vascular": "#bab0ab",
    "Glial": "#b07aa1",
    "Glioma": "#d37295",
}


def save_multi(fig, stem: str):
    fig.savefig(ROOT / f"{stem}.jpg", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
    basis = basis.loc[CELLTYPE_ORDER, ECO_ORDER]

    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    x = np.arange(len(ECO_ORDER))
    bottom = np.zeros(len(ECO_ORDER))

    for ct in CELLTYPE_ORDER:
        vals = basis.loc[ct, ECO_ORDER].to_numpy()
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=CELLTYPE_COLORS[ct],
            edgecolor="white",
            linewidth=0.4,
            width=0.72,
            label=ct,
        )
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([ECO_SHORT[e] for e in ECO_ORDER], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Composition weight", fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.set_title("TLS Ecotype Composition", fontsize=10, pad=8)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.98, 0.5),
        fontsize=7,
        ncol=1,
    )
    fig.tight_layout(rect=[0, 0, 0.84, 1])
    save_multi(fig, "fig_nmf_basis_stacked_bar")
    print("Saved stacked composition barplot.")


if __name__ == "__main__":
    main()
