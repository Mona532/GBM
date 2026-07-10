"""TLS component clustering: group TLS spots into components, cluster by composition"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

SKIP = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
OUT = Path(r"E:/GBM/results"); OUT.mkdir(parents=True, exist_ok=True)
def save(fig,name):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(OUT/f"{name}.{fmt}",bbox_inches="tight",dpi=600)

# Cell types for TLS component characterization
TLS_FEATURES = ["ILC1","ILC2","ILC3","B","CD8_T","Dendritic","NK","Myeloid",
                "AC-gliosis-like","Astrocytes","Vascular-associated","Hypoxic","Proliferative"]

# Collect TLS components across all samples
components = []  # [{sample, size, features: {cell_type: mean_q05}, tls_score_mean}]

for tls_dir, h5_dir in DS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir() or sd.name in SKIP: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"; h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        if "barcode" in tls.columns: tls = tls.set_index("barcode")
        adata = ad.read_h5ad(h5)
        shared = adata.obs_names.intersection(tls.index)
        if len(shared) < 100: continue
        adata = adata[shared]; tls = tls.loc[shared]
        tls_mask = (tls["TLS.region"]=="TLS").values
        n_tls = tls_mask.sum()
        if n_tls < 5: continue

        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        coords = adata.obsm["spatial"]
        tls_idx = np.where(tls_mask)[0]
        tls_coords = coords[tls_idx]

        # Spatial connectivity within TLS spots
        k = min(7, len(tls_idx))
        adj = kneighbors_graph(tls_coords, n_neighbors=k, mode="connectivity", include_self=True)
        n_comp, labels = connected_components(adj, directed=False)

        for cid in range(n_comp):
            comp_mask = labels == cid
            comp_global_idx = tls_idx[comp_mask]
            comp_size = comp_mask.sum()
            if comp_size < 3: continue
            feat = {}
            for ct in TLS_FEATURES:
                if ct in ct_names:
                    feat[ct] = q05[comp_global_idx, ct_names.index(ct)].mean()
            feat["size"] = comp_size
            feat["sample"] = sd.name
            if "TLS.score" in tls.columns:
                feat["tls_score"] = tls.iloc[comp_global_idx]["TLS.score"].mean()
            components.append(feat)

print(f"TLS components: {len(components)}")
df_comp = pd.DataFrame(components)
df_comp.to_csv(OUT / "tls_components.csv", index=False)

# Feature matrix for clustering
feat_cols = TLS_FEATURES + (["tls_score"] if "tls_score" in df_comp.columns else [])
X = df_comp[feat_cols].fillna(0).values
X_scaled = StandardScaler().fit_transform(X)

# PCA for visualization
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Hierarchical clustering
n_clusters = 5
hc = AgglomerativeClustering(n_clusters=n_clusters)
labels = hc.fit_predict(X_scaled)
df_comp["cluster"] = labels

# Per-cluster composition
cluster_means = df_comp.groupby("cluster")[feat_cols].mean()
print(f"\nCluster sizes: {df_comp['cluster'].value_counts().sort_index().to_dict()}")
print("\nPer-cluster mean composition:")
print(cluster_means.to_string())

# Figure: PCA + composition heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

# PCA
colors = plt.cm.tab10.colors[:n_clusters]
for ci in range(n_clusters):
    mask = labels == ci
    ax1.scatter(X_pca[mask,0], X_pca[mask,1], c=[colors[ci]], s=3, alpha=0.5, label=f"C{ci} (n={mask.sum()})")
ax1.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", fontsize=7)
ax1.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", fontsize=7)
ax1.set_title(f"TLS component clusters (n={len(df_comp)})", fontsize=8, fontweight="bold", loc="left")
ax1.legend(fontsize=5.5, loc="upper right")

# Composition heatmap (drop tls_score for z-score)
feat_no_score = [f for f in feat_cols if f != "tls_score"]
comp_mat = cluster_means[feat_no_score].values
comp_z = (comp_mat - comp_mat.mean(axis=1, keepdims=True)) / (comp_mat.std(axis=1, keepdims=True) + 1e-8)
im = ax2.imshow(comp_z.T, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax2.set_xticks(range(n_clusters)); ax2.set_xticklabels([f"C{c}" for c in range(n_clusters)], fontsize=7)
ax2.set_yticks(range(len(feat_no_score))); ax2.set_yticklabels(feat_no_score, fontsize=6)
ax2.set_title("Cluster composition (row z-score)", fontsize=8, fontweight="bold", loc="left")
plt.colorbar(im, ax=ax2, shrink=0.8)

fig.tight_layout(); save(fig, "fig_tls_component_clusters"); plt.close()
print(f"\nFigure: {OUT}/fig_tls_component_clusters.{{svg,pdf,tiff}}")
