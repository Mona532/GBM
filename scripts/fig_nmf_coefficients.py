"""
NMF H matrix: component × ecotype heatmap, sorted by dominant ecotype.
Shows soft assignment — components can belong to multiple ecotypes.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

# Load unit weights (H matrix: component × ecotype)
w = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
eco_cols = ["E1","E2","E3","E4"]
H = w[eco_cols].values

# Sort by dominant ecotype, then by weight within each ecotype
dominant = np.argmax(H, axis=1)
max_weight = H.max(axis=1)
order = np.lexsort((max_weight, dominant))  # stable sort by eco then weight

H_sorted = H[order]
dominant_sorted = dominant[order]

# Color palette
palette = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]
eco_short = ["Lymphocyte", "ILC-enriched", "Myeloid-vascular", "Glial-CD4"]
eco_n = [(dominant_sorted == i).sum() for i in range(4)]

# ====== Figure ======
fig, ax = plt.subplots(figsize=(5, 8))
im = ax.imshow(H_sorted.T, aspect="auto", cmap="YlOrRd", vmin=0, vmax=H_sorted.max() * 0.9)

# Ecotype labels
ax.set_yticks(range(4))
ax.set_yticklabels([f"{eco_short[i]} (n={eco_n[i]})" for i in range(4)], fontsize=8)

# Add boundary lines between ecotypes
boundary = 0
for i in range(4):
    boundary += eco_n[i]
    ax.axhline(y=i+0.5, color="black", linewidth=0.5)
    if i < 3:
        ax.axvline(x=boundary-0.5, color="black", linewidth=0.8, linestyle="--")

# Color bar for ecotype groups at top
for i in range(4):
    start = sum(eco_n[:i])
    end = start + eco_n[i]
    ax.axvspan(start-0.5, end-0.5, ymin=1.02, ymax=1.06, color=palette[i], clip_on=False, alpha=0.8)

ax.set_xticks([])
ax.set_xlabel(f"{len(w)} TLS components", fontsize=8)

cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
cbar.set_label("NMF weight", fontsize=7)

fig.savefig(ROOT / "fig_nmf_coefficients.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {ROOT}/fig_nmf_coefficients.jpg")
print(f"Ecotype sizes: {eco_n}")
