"""
fig_deg_heatmap.py — DEG heatmap for TLS programs

Reads posdev_DEG_heatmap_data.csv (from 03_deg_heatmap.R) and renders a
row-z-scored heatmap of top differentially expressed genes per program.
Only programs with ≥20 significant DEGs (FDR<0.1, logFC>0.5) are shown.
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Nature-style defaults
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.6,
})

root = "E:/GBM/results"
df   = pd.read_csv(f"{root}/posdev_DEG_heatmap_data.csv")

# composition-based names (old CSV columns → new display names)
rename = {
    "Myeloid": "Myeloid",
    "HypoxicVascular": "T-cell",
    "Oligodendrocyte": "Lymphoid-Glial",
    "Inflammatory": "ILC-NK",
    "Neuronal": "TLS-structural",
}
progs_display = ["Myeloid", "T-cell", "Lymphoid-Glial", "ILC-NK", "TLS-structural"]

df["program"] = df["program"].map(rename)
df["program"] = pd.Categorical(df["program"], categories=progs_display, ordered=True)
df = df.sort_values("program")

# row z-score within each gene (use old CSV column names, rename after)
M = df.set_index("gene")[list(rename.keys())]
M.columns = [rename[c] for c in M.columns]
M = M[progs_display]
Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=0).replace(0, 1), axis=0)

# block boundaries for program separators
grp   = df["program"].tolist()
ends  = [grp.count(p) for p in progs_display]
cum   = np.cumsum(ends)
starts = [0] + list(cum[:-1])

fig, ax = plt.subplots(figsize=(4.5, 0.18 * Z.shape[0] + 1.2))
im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")

ax.set_xticks(range(len(progs_display)))
ax.set_xticklabels(progs_display, rotation=30, ha="right", fontsize=7)
ax.set_yticks(range(Z.shape[0]))
ax.set_yticklabels(Z.index, fontsize=5.5)

# horizontal separators and black boxes per program block
for e in cum[:-1]:
    ax.axhline(e - 0.5, color="black", lw=0.8)
for j in range(len(progs_display)):
    ax.add_patch(Rectangle(
        (j - 0.5, starts[j] - 0.5), 1, ends[j],
        fill=False, ec="black", lw=1.4))

ax.tick_params(length=2, width=0.5)
for s in ax.spines.values():
    s.set_linewidth(0.6)

cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label("logCPM (row z-score)", fontsize=6.5)
cb.ax.tick_params(labelsize=6, length=2)
cb.outline.set_linewidth(0.6)

ax.set_title("Top DEGs per program (FDR<0.1, logFC>0.5)",
             fontsize=8, fontweight="bold", pad=6)
fig.savefig(f"{root}/fig_program_topDEG.jpg", dpi=400, bbox_inches="tight")
fig.savefig(f"{root}/fig_program_topDEG.pdf", bbox_inches="tight")
print("Saved fig_program_topDEG.jpg / .pdf")
