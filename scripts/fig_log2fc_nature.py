"""Nature horizontal bar chart — no clipping, clean layout"""
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.5,
    "legend.frameon": False,
})

CAT = {
    "Excit (Glutamate)":  "#E64A19",
    "Inhib (GABA/Gly)":   "#2E7D32",
    "Cholinergic (ACh)":  "#1565C0",
    "DA/NE":              "#7B1FA2",
    "Serotonin (5-HT)":   "#C62828",
}

df = pd.read_csv(r"E:/GBM/results/receptor_tls_log2fc_ilc_samples.csv")
df = df[(df["median_log2FC"] > 1.0) & (df["fdr"] < 0.1) & (df["n_samples"] >= 10)]
df = df.nlargest(10, "median_log2FC").sort_values("median_log2FC")
N = len(df)

fig, ax = plt.subplots(figsize=(6, 4))

# Lollipop stems (thin lines)
for i, (_, r) in enumerate(df.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax.plot([0, r["median_log2FC"]], [i, i], color=c, linewidth=1.5, alpha=0.6,
            solid_capstyle="round", zorder=3)

# Dots
for i, (_, r) in enumerate(df.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax.scatter(r["median_log2FC"], i, s=60, c=c, alpha=0.9, edgecolors="black", linewidths=0.5, zorder=5)

ax.set_yticks(range(N))
ax.set_yticklabels(df["gene"].values, fontsize=8, fontstyle="italic")
ax.set_xlabel("log2(TLS / nonTLS)", fontsize=8, labelpad=6)
ax.set_xlim(0, df["median_log2FC"].max() * 1.2)
ax.invert_yaxis()

# Category legend — inside plot, top-left
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
cats_present = set(df["category"].values)
cat_h = [Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_present]
leg1 = ax.legend(handles=cat_h, fontsize=6, loc="lower right", title="Category", title_fontsize=6.5,
                 handlelength=1, borderpad=0.4, labelspacing=0.3)
ax.add_artist(leg1)

# FDR legend — below first legend
sig_h = [Line2D([0],[0], marker="o", color="w", markerfacecolor="gray", markersize=8, markeredgecolor="black", markeredgewidth=0.5)]
ax.legend(handles=sig_h, fontsize=6, loc="lower right", bbox_to_anchor=(1.0, -0.18),
          title="FDR < 0.1", title_fontsize=6.5, borderpad=0.4)

ax.set_title("Receptor expression: TLS vs non-TLS\n(GBM samples with ILC-high TLS)",
             fontsize=9, fontweight="bold", loc="left", pad=12)

fig.tight_layout(pad=1.5)
fig.savefig(r"E:/GBM/results/fig_receptor_tls_log2fc.svg", bbox_inches="tight")
fig.savefig(r"E:/GBM/results/fig_receptor_tls_log2fc.pdf", bbox_inches="tight")
fig.savefig(r"E:/GBM/results/fig_receptor_tls_log2fc.tiff", dpi=600, bbox_inches="tight")
plt.close()
print(f"Genes: {N} — saved")
