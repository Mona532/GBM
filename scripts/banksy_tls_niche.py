"""
BANKSY spatial niche analysis on all TLS regions (consolidated 13-type reference).
Computes spatial niches from c2l abundance, then quantifies TLS enrichment per niche.
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from sklearn.cluster import KMeans
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

H5AD = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated")
TLS_DIR = Path(r"E:/GBM/results/tls_consolidated")
OUT = Path(r"E:/GBM/results"); OUT.mkdir(parents=True, exist_ok=True)

K = 5; LAMBDA = 0.2; N_NEIGH = 15
ILC_TYPES = ["ILC1","ILC2","ILC3"]

# ====== Step 1: Collect per-sample data ======
data_list, tls_masks, sample_ids = [], [], []
for h5_path in sorted(H5AD.glob("*.h5ad")):
    tls_csv = TLS_DIR / h5_path.stem / "tls_spot_scores_official_relaxed.csv"
    if not tls_csv.exists(): continue
    tls = pd.read_csv(tls_csv)
    if "barcode" in tls.columns: tls = tls.set_index("barcode")
    a = ad.read_h5ad(h5_path)
    shared = a.obs_names.intersection(tls.index)
    if len(shared) < 100: continue
    a = a[shared]; tls = tls.loc[shared]
    tls_mask = (tls["TLS.region"]=="TLS").values
    if tls_mask.sum() < 5: continue

    q05 = a.obsm["c2l_ilc_q05"]
    ct_list = list(a.uns["c2l_ilc_cell_types"])
    coords = a.obsm["spatial"]

    # Build AnnData for BANKSY: features = c2l abundance
    ba = ad.AnnData(X=q05.astype(np.float64),
                     obs=pd.DataFrame({"sample": h5_path.stem}, index=range(q05.shape[0])),
                     obsm={"spatial": coords, "spatial_coords": coords})
    ba.var_names = ct_list
    data_list.append(ba); tls_masks.append(tls_mask); sample_ids.append(h5_path.stem)

print(f"Samples: {len(data_list)}, TLS spots: {sum(m.sum() for m in tls_masks)}")

# ====== Step 2: BANKSY per sample, aggregate ======
niche_tls_counts = np.zeros(K)       # TLS spots per niche (across all samples)
niche_all_counts = np.zeros(K)       # all spots per niche
niche_q05_sum = np.zeros((K, len(ct_list)))  # summed q05 per niche
niche_sample_counts = np.zeros(K)    # number of samples with this niche

for i, (ba, tls_mask, sid) in enumerate(zip(data_list, tls_masks, sample_ids)):
    try:
        bd = initialize_banksy(ba, coord_keys=("spatial","spatial","spatial_coords"),
                               num_neighbours=N_NEIGH, max_m=0,
                               plt_edge_hist=False, plt_nbr_weights=False, plt_theta=False)
        _, bm = generate_banksy_matrix(ba, bd, [LAMBDA], max_m=0)
        X = bm.X.toarray() if hasattr(bm.X, "toarray") else bm.X
        labels = KMeans(n_clusters=K, random_state=42, n_init=10).fit_predict(X)
    except Exception as e:
        print(f"  BANKSY failed for {sid}: {e}")
        continue

    for k in range(K):
        mask_k = labels == k
        niche_tls_counts[k] += (labels[tls_mask] == k).sum()
        niche_all_counts[k] += mask_k.sum()
        niche_q05_sum[k] += ba.X[mask_k].sum(axis=0)
        if mask_k.sum() > 0:
            niche_sample_counts[k] += 1

    if (i+1) % 40 == 0:
        print(f"  {i+1}/{len(data_list)}")

# ====== Step 3: Per-niche metrics ======
niche_comp = niche_q05_sum / (niche_all_counts[:, None] + 1e-8)
global_mean = np.vstack([d.X for d in data_list]).mean(axis=0) + 1e-8
niche_enrich = np.log2(niche_comp / global_mean)

# TLS enrichment per niche
tls_pct = niche_tls_counts / niche_tls_counts.sum() * 100
all_pct = niche_all_counts / niche_all_counts.sum() * 100
total_tls = niche_tls_counts.sum()
total_all = niche_all_counts.sum()
tls_log2fc = np.log2((niche_tls_counts + 1) / (niche_all_counts + 1) * total_all / total_tls)

# ILC features per niche
ilc_frac_high = np.zeros((K, len(ILC_TYPES)))
ilc_means = np.zeros((K, len(ILC_TYPES)))
for k in range(K):
    for j, ct in enumerate(ILC_TYPES):
        idx = ct_list.index(ct) if ct in ct_list else -1
        if idx >= 0:
            ilc_means[k, j] = niche_comp[k, idx]

print(f"\nNiche results (K={K}):")
for k in range(K):
    top_ct = ct_list[np.argmax(niche_enrich[k])]
    print(f"  N{k}: spots={niche_all_counts[k]:.0f}, TLS%={tls_pct[k]:.1f}%, "
          f"TLS_log2FC={tls_log2fc[k]:+.2f}, top={top_ct}, samples={niche_sample_counts[k]:.0f}")

# ====== Figure ======
fig = plt.figure(figsize=(12, 6))

# --- Panel A: Niche composition log2 enrichment ---
ax1 = fig.add_axes([0.05, 0.12, 0.38, 0.80])
vmax = max(abs(niche_enrich.min()), abs(niche_enrich.max()), 0.5)
im1 = ax1.imshow(niche_enrich.T, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax1.set_xticks(range(K))
tls_labels = [f"N{k}\nn={niche_all_counts[k]:.0f}\nTLS {tls_log2fc[k]:+.2f}" for k in range(K)]
ax1.set_xticklabels(tls_labels, fontsize=6)
ax1.set_yticks(range(len(ct_list)))
ax1.set_yticklabels(ct_list, fontsize=6)
ax1.set_title(f"a  BANKSY niche composition (K={K}, lambda={LAMBDA})", fontsize=9, fontweight="bold", loc="left")
plt.colorbar(im1, ax=ax1, shrink=0.7).set_label("log2 enrichment", fontsize=6)

# --- Panel B: TLS enrichment per niche ---
ax2 = fig.add_axes([0.50, 0.12, 0.22, 0.80])
colors = ["#d7191c" if v > 0 else "#2b83ba" for v in tls_log2fc]
bars = ax2.bar(range(K), tls_log2fc, color=colors, alpha=0.7, edgecolor="black", linewidth=0.3)
ax2.axhline(y=0, color="black", linewidth=0.5, linestyle="--")
ax2.set_xticks(range(K))
ax2.set_xticklabels([f"N{k}" for k in range(K)], fontsize=7)
ax2.set_ylabel("TLS enrichment (log2 FC)", fontsize=7)
ax2.set_title("b  TLS enrichment by niche", fontsize=9, fontweight="bold", loc="left")
for bar, v in zip(bars, tls_log2fc):
    ax2.text(bar.get_x() + bar.get_width()/2, v + 0.05*np.sign(v), f"{v:+.2f}", ha="center", fontsize=7)

# --- Panel C: ILC abundance per niche ---
ax3 = fig.add_axes([0.78, 0.12, 0.20, 0.80])
ilc_z = (ilc_means - ilc_means.mean(axis=0)) / (ilc_means.std(axis=0) + 1e-8)
im3 = ax3.imshow(ilc_z.T, aspect="auto", cmap="YlOrRd", vmin=-1, vmax=2)
ax3.set_xticks(range(K))
ax3.set_xticklabels([f"N{k}" for k in range(K)], fontsize=7)
ax3.set_yticks(range(len(ILC_TYPES)))
ax3.set_yticklabels(ILC_TYPES, fontsize=6)
ax3.set_title("c  ILC abundance (z-score)", fontsize=9, fontweight="bold", loc="left")
plt.colorbar(im3, ax=ax3, shrink=0.7).set_label("z-score", fontsize=6)

fig.savefig(OUT / "fig_banksy_tls_niche.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nFigure: {OUT}/fig_banksy_tls_niche.jpg")

# ====== Save data ======
niche_df = pd.DataFrame(niche_enrich, columns=ct_list)
niche_df["niche"] = range(K)
niche_df["tls_log2fc"] = tls_log2fc
niche_df["tls_pct"] = tls_pct
niche_df["n_spots"] = niche_all_counts
niche_df["n_samples"] = niche_sample_counts
niche_df.to_csv(OUT / "banksy_niche_tls.csv", index=False)
print(f"Data: {OUT}/banksy_niche_tls.csv")

# ====== Summary ======
print("\nKey findings:")
for k in range(K):
    top3_idx = np.argsort(-niche_enrich[k])[:3]
    top3 = [f"{ct_list[i]}({niche_enrich[k,i]:+.2f})" for i in top3_idx]
    print(f"  N{k}: TLS_FC={tls_log2fc[k]:+.2f}, top={', '.join(top3)}")
