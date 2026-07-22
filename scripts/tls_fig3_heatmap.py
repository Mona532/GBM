"""Fig3-style marker gene heatmap per NMF ecotype, following Pan-Cancer Atlas."""
import pandas as pd, numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"svg.fonttype":"none","pdf.fonttype":42,
    "font.size":6,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.3,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")
markers = pd.read_csv(r"E:/GBM/Pan-Cancer_Spatial_Atlas_TLS-main/data/fig3_markers.csv")

# Load pseudobulk gene expression per component
meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")
nmf = pd.read_csv(ROOT / "tls_compnmf_rank5_unit_weights.csv")
spot_map = pd.read_csv(ROOT / "tls_component_spot_map.csv")

eco_names = {"E1":"Glial-CD4","E2":"TLS-structural","E3":"Vascular","E4":"Lymphocyte","E5":"Myeloid"}
eco_order = ["Lymphocyte","TLS-structural","Glial-CD4","Vascular","Myeloid"]

# All marker genes
all_genes = markers["gene"].tolist()
# Check against available genes
import anndata as ad
test = ad.read_h5ad("E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated_tls16/AT10-BRA-5-FO-1_1.h5ad")
ge_test = test[:, test.var["feature_types"]=="Gene Expression"]
avail = [g for g in all_genes if g in ge_test.var_names]
missing = set(all_genes) - set(avail)
print(f"Available: {len(avail)}/{len(all_genes)}, missing: {missing}")

# Per-component mean log1p expression per gene
comp_expr = {}
for _, row in meta.iterrows():
    uid = row["unit_id"]; sid = row["sample"]
    h5_path = f"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated_tls16/{sid}.h5ad"
    try:
        a = ad.read_h5ad(h5_path)
        ge = a[:, a.var["feature_types"]=="Gene Expression"]
        spots = spot_map[spot_map["unit_id"]==uid]["barcode"].values
        spots = [s for s in spots if s in a.obs_names]
        if len(spots) < 5: continue
        ge_sub = ge[spots]
        X = ge_sub.X.toarray() if hasattr(ge_sub.X, "toarray") else ge_sub.X
        X_log = np.log1p(X)
        g2i = {g: i for i, g in enumerate(ge_sub.var_names)}
        comp_expr[uid] = {g: np.mean(X_log[:, g2i[g]]) for g in avail if g in g2i}
    except: pass

expr_df = pd.DataFrame(comp_expr).T  # components x genes
expr_df.index.name = "unit_id"
expr_df = expr_df.reset_index()
expr_df = expr_df.merge(nmf[["unit_id","dominant_ecotype"]], on="unit_id", how="left")
expr_df["eco_name"] = expr_df["dominant_ecotype"].map(eco_names)

# Per-ecotype mean, then z-score across ecotypes
eco_means = expr_df.groupby("eco_name")[avail].mean()
eco_means = eco_means.loc[[e for e in eco_order if e in eco_means.index]]
eco_z = StandardScaler().fit_transform(eco_means.T).T  # z-score per gene across ecotypes

# Organize by marker categories
major_types = ["B","T","CAF_Macrophage_DC","Others"]
cat_colors = {"B":"#1B5E20","T":"#E64A19","CAF_Macrophage_DC":"#0D47A1","Others":"#6A1B9A"}
gene_order = []
gene_cats = []
gene_subcats = []
for mt in major_types:
    sub = markers[markers["Major_type"]==mt]
    for ct in sub["celltype"].unique():
        genes = [g for g in sub[sub["celltype"]==ct]["gene"].tolist() if g in avail]
        gene_order.extend(genes)
        gene_cats.extend([mt]*len(genes))
        gene_subcats.extend([ct]*len(genes))

mat = eco_z[:, [avail.index(g) for g in gene_order]]
vmax = max(abs(mat.min()), abs(mat.max()))

fig, ax = plt.subplots(figsize=(20, 4))
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_yticks(range(len(eco_order)))
ax.set_yticklabels(eco_order, fontsize=8)
ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, fontsize=5, rotation=90)

# Category color bar on top
prev = ""
for i, (g, cat) in enumerate(zip(gene_order, gene_cats)):
    if cat != prev:
        ax.axvline(x=i-0.5, color="black", linewidth=0.5)
        prev = cat
    ax.text(i, -1.2, cat[0], fontsize=5, color=cat_colors[cat], ha="center", fontweight="bold")

plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02).set_label("z-score", fontsize=7)
ax.set_title("Marker gene expression by TLS ecotype (Pan-Cancer Atlas Fig3)", fontsize=9, fontweight="bold")
fig.tight_layout()
fig.savefig(ROOT / "fig_tls_fig3_heatmap.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {ROOT}/fig_tls_fig3_heatmap.jpg")
