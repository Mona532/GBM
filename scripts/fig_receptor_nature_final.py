"""Nature-style figure: neurotransmitter receptors in ILC-dominant TLS spots. Single hero panel — only Visium-detectable genes."""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Nature defaults ──
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.5,
    "legend.frameon": False,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

def save_pub(fig, stem, dpi=600):
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.tiff", dpi=dpi, bbox_inches="tight")

# ── Data ──
df = pd.read_csv(r"E:/GBM/results/receptor_gbm_all_per_gene.csv")
df_plot = df[df["log2FC"] > -10].copy()  # only Visium-detectable
df_plot = df_plot.sort_values("log2FC")

# Category palette — restrained, Nature-compatible
cat_colors = {
    "Excit (Glutamate)":    "#c44e52",
    "Inhib (GABA/Gly)":     "#55a868",
    "Cholinergic (ACh)":    "#4c72b0",
    "DA/NE":                "#937860",
    "Serotonin (5-HT)":     "#ccb974",
}

# ── Figure ──
fig, ax = plt.subplots(figsize=(4.2, 3.6))

N = len(df_plot)
y = np.arange(N)

# Reference line
ax.axvline(x=0, color="black", linewidth=0.4, zorder=1)

# Lollipop stems
for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr"] < 0.1
    lw = 1.5 if is_sig else 0.7
    alpha = 0.7 if is_sig else 0.25
    ax.plot([0, r["log2FC"]], [i, i], color=color, linewidth=lw, alpha=alpha,
            solid_capstyle="round", zorder=2)

# Dots
for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr"] < 0.1
    ax.scatter(r["log2FC"], i,
               s=42 if is_sig else 24,
               c=color,
               alpha=0.9 if is_sig else 0.35,
               edgecolors="black" if is_sig else "none",
               linewidths=0.5 if is_sig else 0,
               zorder=4)

# Gene labels (italic, on the bars themselves or right-aligned)
for i, (_, r) in enumerate(df_plot.iterrows()):
    xi = r["log2FC"]
    if xi >= 0:
        ax.text(xi + 0.06, i, r["gene"], ha="left", va="center",
                fontsize=6.5, fontstyle="italic", color="#222222")
    else:
        ax.text(xi - 0.06, i, r["gene"], ha="right", va="center",
                fontsize=6.5, fontstyle="italic", color="#222222")

# Axes
ax.set_yticks([])
ax.set_xlabel("log2(ILC-TLS / non-TLS)", fontsize=7.5, labelpad=6)
ax.set_xlim(-2.2, 2.0)
ax.tick_params(labelsize=6)

# FDR legend (minimal)
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#444444", markersize=7,
           markeredgecolor="black", markeredgewidth=0.5),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#bbbbbb", markersize=5,
           markeredgecolor="none"),
]
ax.legend(legend_elements, ["FDR < 0.1", "not significant"],
          fontsize=6, loc="lower right", handletextpad=0.4, borderpad=0.3,
          labelspacing=0.3)

# Category swatches (compact, bottom-left)
for j, (cat, color) in enumerate(cat_colors.items()):
    ax.text(0.02, 0.98 - j * 0.07, cat, transform=ax.transAxes,
            fontsize=5.5, color=color, va="top", fontweight="bold")

# Title
ax.set_title("Neurotransmitter receptors in ILC-dominant TLS spots",
             fontsize=8, fontweight="bold", loc="left", pad=8)

# Stats line
ax.text(1.0, -0.16, f"n = {df_plot['n_samples'].max()} GBM samples with ILC-dominant TLS  |  13/15 genes FDR < 0.1",
        transform=ax.transAxes, fontsize=5.5, color="#666666", ha="right")

fig.tight_layout(pad=0.5)
save_pub(fig, r"E:/GBM/results/fig_receptor_nature_final")
plt.close()
print("Done: fig_receptor_nature_final.{svg,pdf,tiff}")
