"""BANKSY niche composition + TLS enrichment — save all data"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from sklearn.cluster import KMeans
import warnings; warnings.filterwarnings("ignore")

SKIP = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
OUT = Path(r"E:/GBM/results"); OUT.mkdir(parents=True, exist_ok=True)
K = 5

data_list, tls_masks, all_q05_raw = [], [], []
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
        if tls_mask.sum() < 5: continue
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        coords = adata.obsm["spatial"]
        a = ad.AnnData(X=q05, obs=pd.DataFrame({"sample":[sd.name]*q05.shape[0]}, index=adata.obs_names),
                       obsm={"spatial":coords, "spatial_coords":coords})
        a.var_names = ct_names
        data_list.append(a); tls_masks.append(tls_mask); all_q05_raw.append(q05)

print(f"{len(data_list)} samples")

# BANKSY per sample, aggregate niche data
tls_counts = np.zeros(K); all_counts = np.zeros(K)
niche_q05_sum = np.zeros((K, len(ct_names)))
all_labels = []  # store labels per sample for composition

for i, (a, tls_mask, q05) in enumerate(zip(data_list, tls_masks, all_q05_raw)):
    bd = initialize_banksy(a, coord_keys=("spatial","spatial","spatial_coords"),
                           num_neighbours=15, max_m=0, plt_edge_hist=False,
                           plt_nbr_weights=False, plt_theta=False)
    _, bm = generate_banksy_matrix(a, bd, [0.2], max_m=0)
    X = bm.X; X_dense = X.toarray() if hasattr(X,"toarray") else X
    labels = KMeans(n_clusters=K, random_state=42, n_init=10).fit_predict(X_dense)
    for k in range(K):
        mask_k = labels == k
        tls_counts[k] += (labels[tls_mask]==k).sum()
        all_counts[k] += mask_k.sum()
        niche_q05_sum[k] += q05[mask_k].sum(axis=0)
    all_labels.append(labels)
    if (i+1) % 30 == 0: print(f"  {i+1}/{len(data_list)}")

# Niche composition: mean q05 per niche
niche_comp = niche_q05_sum / (all_counts[:, None] + 1e-8)
global_mean = np.vstack(all_q05_raw).mean(axis=0) + 1e-8
niche_enrich = np.log2(niche_comp / global_mean)

# TLS stats
tls_pct = tls_counts / tls_counts.sum() * 100
all_pct = all_counts / all_counts.sum() * 100
tls_enrich = np.log2((tls_counts + 1) / (all_counts + 1) * all_counts.sum() / tls_counts.sum())

# Save
np.savetxt(OUT / "banksy_niche_composition.csv", niche_enrich, delimiter=",")
pd.DataFrame({"cell_type": ct_names}).to_csv(OUT / "banksy_cell_types.csv", index=False)
pd.DataFrame({"niche": range(K), "tls_pct": tls_pct, "all_pct": all_pct,
              "tls_enrichment": tls_enrich}).to_csv(OUT / "banksy_tls_stats.csv", index=False)
print(f"Saved. Niches={K}")
for ni in range(K):
    top_ct = ct_names[np.argmax(niche_enrich[ni])]
    print(f"  N{ni}: top={top_ct:20s} TLS={tls_pct[ni]:.1f}% enrich={tls_enrich[ni]:+.2f}")
