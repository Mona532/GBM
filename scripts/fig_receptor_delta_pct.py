"""Nature figure: receptor gene positive-rate enrichment in TLS vs non-TLS (per-sample delta_pct)"""
import pandas as pd, numpy as np
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
    for fmt in ["svg", "pdf", "tiff"]:
        fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=dpi)

# ── Data ──
df = pd.read_csv(r"E:/GBM/results/receptor_aggregated_stats.csv")
# Show genes with |delta_pct| > 0.01 — gives ~35 genes, manageable height
THRESH = 0.01
df_plot = df[df["median_delta_pct"].abs() > THRESH].copy()
df_plot = df_plot.sort_values("median_delta_pct")
print(f"Genes shown: {len(df_plot)} (|Δpct| > {THRESH})")

# ── Category colors ──
cat_colors = {
    "Excit (Glutamate)":  "#c44e52",
    "Inhib (GABA/Gly)":   "#55a868",
    "Cholinergic (ACh)":  "#4c72b0",
    "DA/NE":              "#937860",
    "Serotonin (5-HT)":   "#ccb974",
}

# ── Figure dimensions (Nature single column ~89mm = 3.5in) ──
n_genes = len(df_plot)
fig_h = max(3.0, n_genes * 0.22)  # 0.22in per gene
fig, ax = plt.subplots(figsize=(4.5, fig_h))

y = np.arange(n_genes)

# Reference line
ax.axvline(x=0, color="#cccccc", linewidth=0.4, zorder=1, linestyle="--")

# Lollipop stems
for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr_pct"] < 0.1
    lw = 1.5 if is_sig else 0.5
    alpha = 0.7 if is_sig else 0.2
    ax.plot([0, r["median_delta_pct"]], [i, i], color=color, linewidth=lw,
            alpha=alpha, solid_capstyle="round", zorder=2)

# Dots
for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr_pct"] < 0.1
    ax.scatter(r["median_delta_pct"], i, s=42 if is_sig else 20,
               c=color, alpha=0.9 if is_sig else 0.3,
               edgecolors="black" if is_sig else "none",
               linewidths=0.4 if is_sig else 0, zorder=4)

# Gene labels
for i, (_, r) in enumerate(df_plot.iterrows()):
    xi = r["median_delta_pct"]
    if xi >= 0:
        ax.text(xi + 0.006, i, r["gene"], ha="left", va="center",
                fontsize=6.5, fontstyle="italic", color="#222222")
    else:
        ax.text(xi - 0.006, i, r["gene"], ha="right", va="center",
                fontsize=6.5, fontstyle="italic", color="#222222")

# Axes
ax.set_yticks([])
ax.set_xlabel("Δ positive rate (TLS − nonTLS)", fontsize=7.5, labelpad=6)
xmax = max(df_plot["median_delta_pct"].max() * 1.3, 0.1)
xmin = min(df_plot["median_delta_pct"].min() * 1.3, -0.05)
ax.set_xlim(xmin, xmax)
ax.tick_params(labelsize=6)

# FDR legend
ax.legend([
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#444444", markersize=7,
           markeredgecolor="black", markeredgewidth=0.4),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#bbbbbb", markersize=4.5,
           markeredgecolor="none"),
], ["FDR < 0.1", "not significant"], fontsize=6, loc="lower right",
   handletextpad=0.4, borderpad=0.3, labelspacing=0.3)

# Category swatches
for j, (cat, color) in enumerate(cat_colors.items()):
    n = sum(1 for _, r in df_plot.iterrows() if r["category"] == cat)
    if n > 0:
        ax.text(0.02, 0.98 - j * 0.06, f"{cat}  ", transform=ax.transAxes,
                fontsize=5.5, color=color, va="top", fontweight="bold", ha="right")

# Title
ax.set_title("Receptor gene expression: TLS vs non-TLS",
             fontsize=8.5, fontweight="bold", loc="left", pad=8)

# Stats line
n_sig = (df_plot["fdr_pct"] < 0.1).sum()
ax.text(1.0, -0.12,
    f"n = {df_plot['n_samples'].max()} samples  |  {n_sig}/{n_genes} genes FDR < 0.1  |  per-sample Δ positive-rate, Wilcoxon test",
    transform=ax.transAxes, fontsize=5.5, color="#666666", ha="right")

fig.tight_layout(pad=0.8)
save_pub(fig, r"E:/GBM/results/fig_receptor_delta_pct")
plt.close()
print(f"Saved: fig_receptor_delta_pct.{{svg,pdf,tiff}}")
