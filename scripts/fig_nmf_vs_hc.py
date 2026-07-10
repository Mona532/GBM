"""
NMF vs Hierarchical Clustering comparison at K=4.
Shows that NMF resolves distinct ecotypes while HC produces an abundance gradient.
"""
import pandas as pd, numpy as np
from pathlib import Path
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.cluster import AgglomerativeClustering
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

# Load NMF unit weights and composition data
nmf_w = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
comp = pd.read_csv(ROOT / "tls_components_denovo.csv")

# Get cell type proportion columns
prop_cols = [c for c in comp.columns if c.endswith("_prop")]
ct_names = [c.replace("_prop","") for c in prop_cols]

# Merge first, then cluster on the merged subset
comp_nmf = comp.merge(nmf_w[["sample","component_id","dominant_ecotype"]], on=["sample","component_id"], how="inner")
print(f"Merged: {len(comp_nmf)} components (denovo={len(comp)}, nmf={len(nmf_w)})")

# Run HC on the merged subset
prop_cols_avail = [c for c in prop_cols if c in comp_nmf.columns]
X = comp_nmf[prop_cols_avail].to_numpy()
X_clr = np.log(np.clip(X, 1e-8, None))
X_clr = X_clr - X_clr.mean(axis=1, keepdims=True)
X_scaled = RobustScaler().fit_transform(X_clr)
hc_labels = AgglomerativeClustering(n_clusters=4, linkage="ward").fit_predict(X_scaled)
comp_nmf["hc_subtype"] = hc_labels

N_ECO = 4
eco_list = sorted(comp_nmf["dominant_ecotype"].unique())

# ====== Fig 1: Cross-tab heatmap ======
ct_tab = pd.crosstab(comp_nmf["hc_subtype"], comp_nmf["dominant_ecotype"])
# Normalize by row (HC cluster)
ct_norm = ct_tab.div(ct_tab.sum(axis=1), axis=0)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={"width_ratios":[1.2, 1, 1]})

# A: Cross-tab heatmap
im = axes[0].imshow(ct_norm.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
axes[0].set_xticks(range(N_ECO)); axes[0].set_xticklabels(eco_list, fontsize=8)
axes[0].set_yticks(range(4)); axes[0].set_yticklabels([f"HC-{i}" for i in range(4)], fontsize=8)
axes[0].set_xlabel("NMF ecotype", fontsize=8)
axes[0].set_ylabel("HC subtype", fontsize=8)
axes[0].set_title("a  NMF vs HC overlap", fontsize=9, fontweight="bold")
for i in range(4):
    for j in range(N_ECO):
        axes[0].text(j, i, f"{ct_norm.values[i,j]:.2f}", ha="center", va="center", fontsize=7,
                     color="white" if ct_norm.values[i,j] > 0.5 else "black")
plt.colorbar(im, ax=axes[0], shrink=0.7)

# B: HC subtype composition (raw proportion, z-score)
hc_comp = np.zeros((4, len(ct_names)))
for k in range(4):
    hc_comp[k] = comp_nmf[comp_nmf["hc_subtype"]==k][prop_cols].mean().values
hc_z = (hc_comp - hc_comp.mean(axis=0)) / (hc_comp.std(axis=0) + 1e-8)
vmax = max(abs(hc_z.min()), abs(hc_z.max()))
im2 = axes[1].imshow(hc_z.T, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
axes[1].set_yticks(range(len(ct_names))); axes[1].set_yticklabels(ct_names, fontsize=6)
axes[1].set_xticks(range(4))
hc_sizes = [(comp_nmf["hc_subtype"]==k).sum() for k in range(4)]
axes[1].set_xticklabels([f"HC-{k}\nn={hc_sizes[k]}" for k in range(4)], fontsize=7)
axes[1].set_title("b  HC subtypes (K=4, CLR)", fontsize=9, fontweight="bold")
plt.colorbar(im2, ax=axes[1], shrink=0.7).set_label("z-score", fontsize=6)

# C: NMF subtype composition (basis z-score)
basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
basis_z = StandardScaler().fit_transform(basis.values.T).T
vmax3 = max(abs(basis_z.min()), abs(basis_z.max()))
im3 = axes[2].imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax3, vmax=vmax3)
axes[2].set_yticks(range(len(ct_names))); axes[2].set_yticklabels(ct_names, fontsize=6)
axes[2].set_xticks(range(N_ECO)); axes[2].set_xticklabels(eco_list, fontsize=8)
axes[2].set_title("c  NMF ecotypes (basis)", fontsize=9, fontweight="bold")
plt.colorbar(im3, ax=axes[2], shrink=0.7).set_label("z-score", fontsize=6)

fig.tight_layout()
fig.savefig(ROOT / "fig_nmf_vs_hc.jpg", dpi=300, bbox_inches="tight")
plt.close()

# Print summary
print("NMF vs HC comparison:")
print(f"\nHC subtype sizes: {hc_sizes}")
print(f"NMF ecotype sizes: {[(comp_nmf['dominant_ecotype']==e).sum() for e in eco_list]}")

# How mixed are HC subtypes vs NMF?
hc_purity = ct_norm.max(axis=1).mean()
nmf_col = pd.crosstab(comp_nmf["hc_subtype"], comp_nmf["dominant_ecotype"])
nmf_col_norm = nmf_col.div(nmf_col.sum(axis=0), axis=0)
nmf_purity = nmf_col_norm.max(axis=0).mean()
print(f"HC purity (mean max row): {hc_purity:.2f}")
print(f"NMF purity (mean max col): {nmf_purity:.2f}")
print(f"Saved: {ROOT}/fig_nmf_vs_hc.jpg")
