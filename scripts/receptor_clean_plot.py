"""Clean figure: only genes with detectable expression in both groups"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT_DIR = Path(r"E:\GBM\results")
df = pd.read_csv(OUT_DIR / "receptor_per_gene_ilc.csv")

# Filter: remove -inf (zero in ILC-TLS), keep expressed in both groups
df_expressed = df[df["log2FC"] > -10].copy()  # keep real FC values
df_expressed = df_expressed.sort_values("log2FC", ascending=False)

rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excitatory (Glutamate)", "Inhibitory (GABA/Gly)", "Cholinergic (ACh)", "DA/NE", "Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna():
        gene_cat[g] = cat_order[idx]

cat_colors = {"Excitatory (Glutamate)": "#E64A19", "Inhibitory (GABA/Gly)": "#1B5E20",
              "Cholinergic (ACh)": "#0D47A1", "DA/NE": "#6A1B9A", "Serotonin (5-HT)": "#BF360C"}

fig, ax = plt.subplots(figsize=(14, 8))
y = np.arange(len(df_expressed))
x = df_expressed["log2FC"].values
err_low = np.maximum(0, x - 0.2)  # placeholder
err_high = x + 0.2

colors = [cat_colors.get(gene_cat.get(g, ""), "gray") for g in df_expressed["gene"]]
markers = ["s" if fdr < 0.05 else "o" for fdr in df_expressed["fdr"]]

# Bar plot
bars = ax.barh(y, x, height=0.7, color=colors, alpha=0.35, edgecolor="none")
# Significant markers
for i, (_, r) in enumerate(df_expressed.iterrows()):
    if r["fdr"] < 0.1:
        ax.scatter(r["log2FC"] + 0.02, i, c=colors[i], s=60, edgecolors="black", linewidths=0.8, zorder=5)

ax.set_yticks(y)
ax.set_yticklabels(df_expressed["gene"].values, fontsize=9)
ax.axvline(x=0, color="black", linewidth=1)
ax.set_xlabel("log2(ILC-TLS / non-TLS)", fontsize=13)
ax.set_title("Detectable neurotransmitter receptors: ILC-dominant TLS vs non-TLS\n(FDR-significant marked with filled dots)", fontsize=14)

# Legend
patches = [mpatches.Patch(color=cat_colors[cn], alpha=0.35, label=cn) for cn in cat_order]
ax.legend(handles=patches, fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_detected_bar.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved to {OUT_DIR}/receptor_detected_bar.png")
print(f"\nExpressed in both groups: {len(df_expressed)} genes")
print(f"FDR<0.1: {(df_expressed['fdr']<0.1).sum()} genes")
