"""Unified-colour Nature figure: receptor enrichment in TLS vs non-TLS (positive-rate delta)"""
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Unified category palette ──
CAT_COLORS = {
    "Excit (Glutamate)":  "#E64A19",
    "Inhib (GABA/Gly)":   "#2E7D32",
    "Cholinergic (ACh)":  "#1565C0",
    "DA/NE":              "#7B1FA2",
    "Serotonin (5-HT)":   "#C62828",
}

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
})

def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]:
        fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

# ── Data ──
df = pd.read_csv(r"E:/GBM/results/receptor_aggregated_stats.csv")
# Only enriched (positive delta_pct), FDR<0.1, and meaningful magnitude
df_pos = df[(df["median_delta_pct"] > 0.005) & (df["fdr_pct"] < 0.1)].copy()
df_pos = df_pos.sort_values("median_delta_pct")
print(f"Genes: {len(df_pos)} enriched (Δpct>0.005, FDR<0.1)")

# ── Figure ──
n = len(df_pos)
fig, ax = plt.subplots(figsize=(4.5, max(3.5, n * 0.22)))

y = np.arange(n)
ax.axvline(x=0, color="black", linewidth=0.4, zorder=1)

for i, (_, r) in enumerate(df_pos.iterrows()):
    color = CAT_COLORS.get(r["category"], "#888888")
    ax.plot([0, r["median_delta_pct"]], [i, i], color=color, linewidth=1.5,
            alpha=0.7, solid_capstyle="round", zorder=2)
    ax.scatter(r["median_delta_pct"], i, s=42, c=color, alpha=0.9,
               edgecolors="black", linewidths=0.4, zorder=4)

for i, (_, r) in enumerate(df_pos.iterrows()):
    xi = r["median_delta_pct"]
    ax.text(xi + 0.005, i, r["gene"], ha="left", va="center",
            fontsize=6.5, fontstyle="italic", color="#222222")

ax.set_yticks([])
ax.set_xlabel("Δ positive rate (TLS − nonTLS)", fontsize=7.5, labelpad=6)
ax.set_xlim(-0.02, df_pos["median_delta_pct"].max() * 1.35)
ax.tick_params(labelsize=6)

# Category legend
for j, (cat, color) in enumerate(CAT_COLORS.items()):
    n_cat = sum(1 for _, r in df_pos.iterrows() if r["category"] == cat)
    if n_cat > 0:
        ax.text(0.98, 0.98 - j * 0.06, f"{cat} ({n_cat})", transform=ax.transAxes,
                fontsize=5.5, color=color, va="top", ha="right", fontweight="bold")

# Title
ax.set_title("Receptor gene expression: TLS vs non-TLS",
             fontsize=8.5, fontweight="bold", loc="left", pad=8)
ax.text(1.0, -0.12,
    f"n = {df_pos['n_samples'].max()} GBM samples  |  all genes FDR < 0.1  |  per-sample Δ positive-rate, Wilcoxon",
    transform=ax.transAxes, fontsize=5.5, color="#666666", ha="right")

fig.tight_layout(pad=0.8)
save_pub(fig, r"E:/GBM/results/fig_receptor_delta_pct")
plt.close()
print(f"Saved: fig_receptor_delta_pct.{{svg,pdf,tiff}}")
