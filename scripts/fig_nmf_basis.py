"""
Single clean NMF basis heatmap — minimal publication style.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_annotated_summary.csv")

ct_list = list(basis.index)
N_ECO = 4

# z-score across ecotypes (row-wise)
basis_z = StandardScaler().fit_transform(basis.values.T).T
vmax = max(abs(basis_z.min()), abs(basis_z.max()))

eco_names = ["Lymphocyte\nTLS", "ILC-enriched\nTLS", "Myeloid-vascular\nTLS", "Glial-CD4\nTLS"]
eco_n = [summary.loc[i, "n_units"] for i in range(N_ECO)]
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]

fig, ax = plt.subplots(figsize=(5, 4.5))
im = ax.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

ax.set_yticks(range(len(ct_list)))
ax.set_yticklabels(ct_list, fontsize=7)

# Ecotype labels with colored text
ax.set_xticks(range(N_ECO))
for i in range(N_ECO):
    ax.text(i, len(ct_list) + 0.3, eco_names[i], ha="center", va="top", fontsize=7,
            fontweight="bold", color=colors[i])
    ax.text(i, len(ct_list) + 1.9, f"n={eco_n[i]}", ha="center", va="top", fontsize=6.5, color="grey")
ax.set_xticks([])
ax.set_xlim(-0.5, N_ECO - 0.5)
ax.set_ylim(len(ct_list) + 2.5, -0.5)

cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("z-score", fontsize=7)

fig.savefig(ROOT / "fig_nmf_basis.jpg", dpi=300, bbox_inches="tight")
fig.savefig(ROOT / "fig_nmf_basis.pdf", bbox_inches="tight")
plt.close()
print(f"Saved: {ROOT}/fig_nmf_basis.jpg")
