"""Figures: CellCharter niches, null receptor result, TLS fraction per niche"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
import cellcharter as cc, squidpy as sq
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

OUT = Path(r"E:/GBM/results")
CAT = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0",
       "DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}

# ==== Fig 1: Receptor null result ====
df = pd.read_csv(OUT / "receptor_cellcharter_final.csv")
df_plot = df.nlargest(20, "median_delta").sort_values("median_delta", ascending=True)

fig, ax = plt.subplots(figsize=(5, 5))
y = np.arange(len(df_plot))
for i, (_, r) in enumerate(df_plot.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.4, solid_capstyle="round", zorder=3)
    ax.scatter(r["median_delta"], i, s=50, c=c, alpha=0.7, edgecolors="#999", linewidths=0.3, zorder=5)
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot["gene"].values, fontsize=7, fontstyle="italic")
ax.axvline(x=0, color="black", linewidth=0.4)
ax.set_xlabel("delta mean log-expr (ILC-enriched - other TLS)", fontsize=8)
ax.set_title("Receptor expression: niche-stratified comparison\n(CellCharter 13 niches, n=57 samples)", fontsize=9, fontweight="bold", loc="left", pad=10)
ax.text(0.98, 0.02, "All FDR > 0.1 - no enrichment detected", transform=ax.transAxes, fontsize=7, color="#c44e52", ha="right", fontweight="bold")
cats_in = set(df_plot["category"].values)
ax.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in],
          fontsize=6, loc="lower right", borderpad=0.3)
fig.tight_layout(pad=0.8)
save_pub(fig, OUT / "fig_null_receptor_niche")
plt.close()
print("Fig 1: receptor null result")

# ==== Fig 2: Niche spatial map + composition on demo sample ====
h5 = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc/AT10-BRA-5-FO-1_2.h5ad")
tls_csv = Path(r"E:/GBM/results/tls_official_cut01/AT10-BRA-5-FO-1_2/tls_spot_scores_official_relaxed.csv")
adata = ad.read_h5ad(h5)
tls = pd.read_csv(tls_csv)
if "barcode" in tls.columns: tls = tls.set_index("barcode")
shared = adata.obs_names.intersection(tls.index)
adata = adata[shared]; tls = tls.loc[shared]
q05 = adata.obsm["c2l_ilc_q05"]
if hasattr(q05,"values"): q05 = q05.values
ct_names = list(adata.uns["c2l_ilc_cell_types"])
coords = adata.obsm["spatial"]
tls_mask = (tls["TLS.region"]=="TLS").values

# CellCharter on demo
obs_df = pd.DataFrame(index=adata.obs_names); obs_df["sample"] = "demo"
a = ad.AnnData(X=q05.astype(np.float32), obs=obs_df, obsm={"spatial": coords})
a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(a, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
autok = cc.tl.ClusterAutoK(n_clusters=(5, 15), max_runs=3)
autok.fit(a, use_rep="X_cellcharter")
niche = autok.predict(a, use_rep="X_cellcharter")
niche_ids = sorted(set(niche))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))

# Spatial map
cmap = ListedColormap(plt.cm.tab10.colors[:len(niche_ids)])
ax1.scatter(coords[:,0], coords[:,1], c=niche.astype(int), s=0.5, cmap=cmap, rasterized=True)
ax1.set_title("AT10-BRA-5-FO-1_2\nSpatial niches (CellCharter)", fontsize=8, fontweight="bold")
ax1.axis("off")

# Niche composition
niche_means = np.array([[q05[niche==n, i].mean() for i in range(len(ct_names))] for n in niche_ids])
im = ax2.imshow(niche_means.T, aspect="auto", cmap="YlOrRd")
ax2.set_xticks(range(len(niche_ids)))
ax2.set_xticklabels([f"N{n}" for n in niche_ids], fontsize=6)
ax2.set_yticks(range(len(ct_names)))
ax2.set_yticklabels(ct_names, fontsize=5)
ax2.set_title("Mean c2l abundance per niche", fontsize=8, fontweight="bold", loc="left")
plt.colorbar(im, ax=ax2, shrink=0.7)
fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_niche_composition_demo")
plt.close()
print("Fig 2: niche composition on demo sample")

# ==== Fig 3: Per-niche TLS fraction across several samples ====
DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
all_niche_tls = {}
count = 0
for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv_d = sd / "tls_spot_scores_official_relaxed.csv"
        h5_d = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv_d.exists() or not h5_d.exists(): continue
        tls_d = pd.read_csv(tls_csv_d)
        if "barcode" in tls_d.columns: tls_d = tls_d.set_index("barcode")
        adata_d = ad.read_h5ad(h5_d)
        shared_d = adata_d.obs_names.intersection(tls_d.index)
        if len(shared_d) < 100: continue
        adata_d = adata_d[shared_d]; tls_d = tls_d.loc[shared_d]
        tls_mask_d = (tls_d["TLS.region"]=="TLS").values
        if tls_mask_d.sum() < 10: continue
        q05_d = adata_d.obsm["c2l_ilc_q05"]
        if hasattr(q05_d,"values"): q05_d = q05_d.values
        ct_d = list(adata_d.uns["c2l_ilc_cell_types"])
        coords_d = adata_d.obsm["spatial"]
        obs_d = pd.DataFrame(index=adata_d.obs_names); obs_d["sample"] = sd.name
        a_d = ad.AnnData(X=q05_d.astype(np.float32), obs=obs_d, obsm={"spatial": coords_d})
        a_d.var_names = ct_d; a_d.obs["sample"] = a_d.obs["sample"].astype("category")
        try:
            sq.gr.spatial_neighbors(a_d, coord_type="generic", delaunay=True)
            cc.gr.aggregate_neighbors(a_d, n_layers=1, use_rep=None, out_key="X_cc", sample_key="sample")
            autok_d = cc.tl.ClusterAutoK(n_clusters=(5, 15), max_runs=3)
            autok_d.fit(a_d, use_rep="X_cc")
            niches_d = autok_d.predict(a_d, use_rep="X_cc")
        except:
            continue
        for n in set(niches_d):
            all_niche_tls.setdefault(n, {"total": 0, "tls": 0})
            all_niche_tls[n]["total"] += (niches_d == n).sum()
            all_niche_tls[n]["tls"] += (niches_d == n).sum()
        count += 1
        if count >= 15: break
    if count >= 15: break

niche_df = pd.DataFrame(all_niche_tls).T
niche_df["pct_tls"] = niche_df["tls"] / niche_df["total"] * 100
niche_df = niche_df.sort_index()

fig, ax = plt.subplots(figsize=(4, 2.5))
x = np.arange(len(niche_df))
ax.bar(x, niche_df["pct_tls"].values, color="#4c72b0", alpha=0.6, edgecolor="#4c72b0", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels([f"N{n}" for n in niche_df.index], fontsize=7)
ax.set_ylabel("% TLS spots", fontsize=7)
ax.set_title("TLS fraction per niche (15 samples)", fontsize=8, fontweight="bold", loc="left")
mean_pct = niche_df["pct_tls"].mean()
ax.axhline(y=mean_pct, color="#c44e52", linestyle="--", linewidth=0.5, label=f"mean={mean_pct:.1f}%")
ax.legend(fontsize=6)
fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_niche_tls_fraction")
plt.close()
print("Fig 3: TLS fraction per niche")
print("All figures saved to:", OUT)
