"""TLS component clustering v2: spatial connectivity + ILC-aware features + robust scaling"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import RobustScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

H5AD = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated")
TLS_DIR = Path(r"E:/GBM/results/tls_consolidated")
OUT = Path(r"E:/GBM/results")

# ====== Parameters ======
K_NEIGH = 7           # spatial neighbors (6 + self)
MIN_SPOTS = 5         # minimum TLS component size
ILC_TYPES = ["ILC1","ILC2","ILC3"]
ILC_THRESH = {"ILC1":1.087,"ILC2":1.077,"ILC3":1.098}  # global TLS P75

# ====== Step 1: Collect TLS components ======
components = []  # [{sample, size, features_dict}]
n_samples = 0

for h5_path in sorted(H5AD.glob("*.h5ad")):
    tls_csv = TLS_DIR / h5_path.stem / "tls_spot_scores_official_relaxed.csv"
    if not tls_csv.exists(): continue
    tls = pd.read_csv(tls_csv)
    if "barcode" in tls.columns: tls = tls.set_index("barcode")
    adata = ad.read_h5ad(h5_path)
    shared = adata.obs_names.intersection(tls.index)
    if len(shared) < 100: continue
    adata = adata[shared]; tls = tls.loc[shared]
    tls_mask = (tls["TLS.region"]=="TLS").values
    if tls_mask.sum() < 5: continue
    n_samples += 1

    q05 = adata.obsm["c2l_ilc_q05"]
    ct_list = list(adata.uns["c2l_ilc_cell_types"])
    coords = adata.obsm["spatial"]
    tls_idx = np.where(tls_mask)[0]
    tls_coords = coords[tls_idx]

    # Spatial connectivity (k-NN on TLS spots only)
    k = min(K_NEIGH, len(tls_idx))
    adj = kneighbors_graph(tls_coords, n_neighbors=k, mode="connectivity", include_self=True)
    n_comp, labels = connected_components(adj, directed=False)

    for cid in range(n_comp):
        mask = labels == cid
        if mask.sum() < MIN_SPOTS: continue
        comp_idx = tls_idx[mask]

        # === Features ===
        feat = {}
        # 1. Mean q05 per cell type
        for ct in ct_list:
            feat[f"{ct}_mean"] = q05[comp_idx, ct_list.index(ct)].mean()

        # 2. ILC-specific features
        total_ilc = 0
        for ct in ILC_TYPES:
            vals = q05[comp_idx, ct_list.index(ct)]
            feat[f"{ct}_mean"] = vals.mean()
            feat[f"{ct}_P90"] = np.percentile(vals, 90)
            feat[f"{ct}_frac_high"] = (vals >= ILC_THRESH[ct]).mean()
            total_ilc += vals
        ilc_total_mean = total_ilc.mean()

        # 3. ILC ratio (ILC / total)
        total_all = sum(q05[comp_idx, ct_list.index(ct)].mean() for ct in ct_list) + 1e-8
        for ct in ILC_TYPES:
            feat[f"{ct}_ratio"] = feat[f"{ct}_mean"] / total_all

        # 4. TLS metadata
        feat["tls_size"] = mask.sum()
        if "TLS.score" in tls.columns:
            feat["tls_score"] = tls.iloc[comp_idx]["TLS.score"].mean()

        feat["sample"] = h5_path.stem
        components.append(feat)

print(f"TLS components: {len(components)} from {n_samples} samples")

# ====== Step 2: Feature matrix + robust scaling ======
df = pd.DataFrame(components)
# Select ILC-aware features for clustering
cluster_cols = (
    [f"{ct}_mean" for ct in ct_list] +
    [f"{ct}_P90" for ct in ILC_TYPES] +
    [f"{ct}_frac_high" for ct in ILC_TYPES] +
    [f"{ct}_ratio" for ct in ILC_TYPES] +
    ["tls_size"]
)
# Only columns that exist
cluster_cols = [c for c in cluster_cols if c in df.columns]
X = df[cluster_cols].fillna(0).astype(float).values

# log1p transform
X_log = np.log1p(X)

# Robust scaling (median + IQR, immune to outliers)
scaler = RobustScaler()
X_scaled = scaler.fit_transform(X_log)

# ====== Step 3: Test K=3-8 ======
print("\nSilhouette scores:")
best_k, best_score = 5, -1
for k in range(3, 9):
    labels = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X_scaled)
    score = silhouette_score(X_scaled, labels)
    print(f"  K={k}: {score:.3f}")
    if score > best_score:
        best_score = score; best_k = k

print(f"\nBest K = {best_k}")

# ====== Step 4: Final clustering ======
labels = KMeans(n_clusters=best_k, random_state=42, n_init=10).fit_predict(X_scaled)
df["cluster"] = labels

# Per-cluster mean composition (log2 enrichment over global)
global_means = {}
for ct in ct_list:
    global_means[ct] = np.mean([c[f"{ct}_mean"] for c in components])
niche_means = np.zeros((best_k, len(ct_list)))
for k in range(best_k):
    mask = labels == k
    for i, ct in enumerate(ct_list):
        niche_means[k, i] = np.log2(df.loc[mask, f"{ct}_mean"].mean() / (global_means[ct] + 1e-8))

# ILC features per cluster
ilc_feats = ["ILC1_frac_high","ILC2_frac_high","ILC3_frac_high","ILC1_ratio","ILC2_ratio","ILC3_ratio","tls_size"]
ilc_feats = [f for f in ilc_feats if f in df.columns]
ilc_by_cluster = df.groupby("cluster")[ilc_feats].mean()

# ====== Figure ======
fig = plt.figure(figsize=(14, 5))

# Composition heatmap
ax1 = fig.add_axes([0.05, 0.15, 0.45, 0.78])
vmax = max(abs(niche_means.min()), abs(niche_means.max()))
im = ax1.imshow(niche_means.T, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax1.set_xticks(range(best_k))
ax1.set_xticklabels([f"C{n}\n({(labels==n).sum()})" for n in range(best_k)], fontsize=7)
ax1.set_yticks(range(len(ct_list))); ax1.set_yticklabels(ct_list, fontsize=6)
ax1.set_title(f"TLS component subtypes (log2 enrichment, {len(components)} components)", fontsize=8, fontweight="bold", loc="left")
plt.colorbar(im, ax=ax1, shrink=0.7)

# ILC features
ax2 = fig.add_axes([0.55, 0.15, 0.42, 0.78])
ilc_mat = ilc_by_cluster.values
ilc_z = (ilc_mat - ilc_mat.mean(axis=0)) / (ilc_mat.std(axis=0) + 1e-8)
im2 = ax2.imshow(ilc_z.T, aspect="auto", cmap="YlOrRd", vmin=-1, vmax=2)
ax2.set_xticks(range(best_k))
ax2.set_xticklabels([f"C{n}" for n in range(best_k)], fontsize=7)
ax2.set_yticks(range(len(ilc_feats))); ax2.set_yticklabels(ilc_feats, fontsize=6)
ax2.set_title("ILC features per cluster (z-score)", fontsize=8, fontweight="bold", loc="left")
plt.colorbar(im2, ax=ax2, shrink=0.7)

fig.savefig(OUT / "fig_tls_component_subtypes_v2.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nFigure: {OUT}/fig_tls_component_subtypes_v2.jpg")
df.to_csv(OUT / "tls_components_v2.csv", index=False)
print(f"Data: {OUT}/tls_components_v2.csv")
print("\nPer-cluster ILC features:")
print(ilc_by_cluster.round(3).to_string())
