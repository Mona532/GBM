"""Nature figure: TLS vs non-TLS log2FC for receptor genes (unified colours)"""
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CAT_COLORS = {
    "Excit (Glutamate)":  "#E64A19",
    "Inhib (GABA/Gly)":   "#2E7D32",
    "Cholinergic (ACh)":  "#1565C0",
    "DA/NE":              "#7B1FA2",
    "Serotonin (5-HT)":   "#C62828",
}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.5,
    "legend.frameon": False,
})

def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]:
        fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

# ── Data ──
df = pd.read_csv(r"E:/GBM/results/receptor_tls_log2fc_ilc_samples.csv")
# Only enriched (positive log2FC), FDR<0.1, enough samples for reliability
df_pos = df[(df["median_log2FC"] > 1.0) & (df["fdr"] < 0.1) & (df["n_samples"] >= 10)].copy()
df_pos = df_pos.nlargest(10, "median_log2FC").sort_values("median_log2FC")
df_pos = df_pos.sort_values("median_log2FC")
print(f"Genes: {len(df_pos)} (top 10, log2FC>1.0, FDR<0.1, n>=10)")

# ── Figure ──
n = len(df_pos)
fig_h = max(3.0, n * 0.32)  # more space per gene
fig, ax = plt.subplots(figsize=(5.0, fig_h))
y = np.arange(n)

ax.axvline(x=0, color="black", linewidth=0.4, zorder=1)

for i, (_, r) in enumerate(df_pos.iterrows()):
    color = CAT_COLORS.get(r["category"], "#888888")
    ax.plot([0, r["median_log2FC"]], [i, i], color=color, linewidth=1.5, alpha=0.7,
            solid_capstyle="round", zorder=2)
    ax.scatter(r["median_log2FC"], i, s=42, c=color, alpha=0.9,
               edgecolors="black", linewidths=0.4, zorder=4)

for i, (_, r) in enumerate(df_pos.iterrows()):
    xi = r["median_log2FC"]
    ax.text(xi + 0.08, i, r["gene"], ha="left", va="center",
            fontsize=6.5, fontstyle="italic", color="#222222")

ax.set_yticks([])
ax.set_xlabel("log2(TLS / nonTLS)", fontsize=7.5, labelpad=6)
xmax = df_pos["median_log2FC"].max() * 1.35  # more breathing room
ax.set_xlim(-0.3, xmax)
ax.tick_params(labelsize=6)

for j, (cat, color) in enumerate(CAT_COLORS.items()):
    n_cat = sum(1 for _, r in df_pos.iterrows() if r["category"] == cat)
    if n_cat > 0:
        ax.text(0.02, 0.98 - j*0.08, f"{cat} ({n_cat})", transform=ax.transAxes,
                fontsize=5.5, color=color, va="top", ha="left", fontweight="bold")

ax.set_title("Receptor expression: TLS vs non-TLS\n(only samples with ILC-high TLS spots)",
             fontsize=8.5, fontweight="bold", loc="left", pad=10)
ax.text(1.0, -0.18,
    f"n = {df_pos['n_samples'].max()} GBM samples with ILC-high TLS  |  per-sample log2FC median, Wilcoxon + FDR",
    transform=ax.transAxes, fontsize=5.5, color="#666666", ha="right")

fig.tight_layout(pad=1.0)
save_pub(fig, r"E:/GBM/results/fig_receptor_tls_log2fc")
plt.close()
print(f"Saved: fig_receptor_tls_log2fc.{{svg,pdf,tiff}}")
