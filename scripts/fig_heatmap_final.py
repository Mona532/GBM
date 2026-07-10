"""
Heatmap-driven final figure: niche composition, ILC-rich TLS distribution, receptor/ligand detection
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
from sklearn.preprocessing import StandardScaler
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})
def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

CAT_COLORS = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
RECEPTOR_GENES = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): RECEPTOR_GENES[g] = list(CAT_COLORS.keys())[idx]
LIGANDS = ["NMU","VIP","ADM","CALCA","CALCB","NPY","TAC1","TAC3","NTS","CCK","GRP","GAL","PENK","PDYN","PNOC","CRH","UCN","AGRP","POMC","MCH"]
ALL_GENES = sorted(set(list(RECEPTOR_GENES.keys()) + LIGANDS))

ILC_TYPES = ["ILC1","ILC2","ILC3"]
GLOBAL_THRESH = {"ILC1":1.034,"ILC2":1.0,"ILC3":1.035}
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results")
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# ====== Step 1: Collect data, run CellCharter, compute per-gene metrics ======
adata_list, sample_data = [], []
for tls_dir, h5_dir in DS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        if "barcode" in tls.columns: tls = tls.set_index("barcode")
        adata = ad.read_h5ad(h5)
        shared = adata.obs_names.intersection(tls.index)
        if len(shared) < 100: continue
        adata = adata[shared]; tls = tls.loc[shared]
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 10: continue
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx_map = {c: ct_names.index(c) for c in ILC_TYPES}
        coords = adata.obsm["spatial"]
        q05_z = (q05 - q05.mean(axis=0)) / (q05.std(axis=0) + 1e-8)
        obs_df = pd.DataFrame(index=adata.obs_names); obs_df["sample"] = sd.name
        a = ad.AnnData(X=q05_z.astype(np.float32), obs=obs_df, obsm={"spatial": coords})
        a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
        adata_list.append(a)
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        vn_all = ge.var_names.values
        genes_present = [g for g in ALL_GENES if g in vn_all]
        gene_idx = [list(vn_all).index(g) for g in genes_present]
        rx_raw = ge[:, gene_idx].X.toarray() if hasattr(ge[:, gene_idx].X,"toarray") else np.asarray(ge[:, gene_idx].X)
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))
        sample_data.append((sd.name, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords))

adata_all = ad.concat(adata_list, join="inner")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adata_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
autok = cc.tl.ClusterAutoK(n_clusters=(5, 20), max_runs=5)
autok.fit(adata_all, use_rep="X_cellcharter")
niche_all = autok.predict(adata_all, use_rep="X_cellcharter")
K = len(set(niche_all))
print(f"Global niches: {K}")

# Niche composition: mean c2l per niche
offset = 0; niche_composition = np.zeros((K, len(ct_names)))
for sid, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords in sample_data:
    n = len(tls_mask); niches = niche_all[offset:offset+n]; offset += n
    for ni in range(K):
        mask = niches == ni
        if mask.sum() > 0:
            niche_composition[ni] += q05_z[mask].sum(axis=0)
ns = np.bincount(niche_all.astype(int), minlength=K) + 1e-8
niche_composition = niche_composition / ns[:, None]

# Per-gene metrics (sample-level)
gene_sample = {}; niche_dist_samples = {}
offset = 0
for sid, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords in sample_data:
    n = len(tls_mask); niches = niche_all[offset:offset+n]; offset += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
    ilc_rich = ilc_high & tls_mask; other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue
    niche_dist_samples[sid] = [int((niches == ni).sum()) for ni in range(K)]
    # ILC-rich distribution across niches
    ilc_niche_counts = [int((niches == ni).sum()) for ni in range(K)]
    ie = expr_norm[ilc_rich]; oe = expr_norm[other_tls]
    for gi, g in enumerate(genes_present):
        delta = ie[:, gi].mean() - oe[:, gi].mean()
        pct_ilc = (rx_raw[ilc_rich, gi] > 0).mean(); pct_other = (rx_raw[other_tls, gi] > 0).mean()
        gene_sample.setdefault(g, []).append({"delta":delta,"pct_ilc":pct_ilc,"pct_other":pct_other})

# Aggregate
rows = []
for g, svals in gene_sample.items():
    n_s = len(svals); deltas = np.array([v["delta"] for v in svals])
    pct_i = np.array([v["pct_ilc"] for v in svals]); pct_o = np.array([v["pct_other"] for v in svals])
    _, p_delta = wilcoxon(deltas, alternative="two-sided")
    rows.append({"gene":g,"n":n_s,"median_delta":np.median(deltas),
        "pct_ILC":np.median(pct_i),"pct_other":np.median(pct_o),
        "prop_positive":(deltas>0).mean(),"p_delta":p_delta,
        "type":"receptor" if g in RECEPTOR_GENES else "ligand",
        "category":RECEPTOR_GENES.get(g,"Ligand")})

df_all = pd.DataFrame(rows)
df_all = df_all[df_all["p_delta"].notna() & (df_all["p_delta"]>0) & (df_all["p_delta"]<1)]
df_all["fdr"] = false_discovery_control(df_all["p_delta"].values)
df_all["candidate"] = ((df_all["n"]>=15)&(df_all["pct_ILC"]>=0.03)&(df_all["median_delta"]>0)&(df_all["prop_positive"]>=0.55))
df_all.to_csv(OUT / "receptor_ligand_heatmap_data.csv", index=False)

# ====== Figure: 3 heatmaps + 1 spatial ======
fig = plt.figure(figsize=(12, 10))

# Panel A: Niche composition heatmap (top-left)
ax_a = fig.add_axes([0.05, 0.55, 0.40, 0.40])
im_a = ax_a.imshow(niche_composition, aspect="auto", cmap="RdBu_r", vmin=-0.3, vmax=0.3)
ax_a.set_xticks(range(len(ct_names))); ax_a.set_xticklabels(ct_names, fontsize=4.5, rotation=90)
ax_a.set_yticks(range(K)); ax_a.set_yticklabels([f"N{n}" for n in range(K)], fontsize=5.5)
ax_a.set_title(f"a  Niche composition (z-scored c2l, {K} niches)", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_a, ax=ax_a, shrink=0.7)

# Panel B: ILC-rich TLS distribution (top-right)
ax_b = fig.add_axes([0.52, 0.55, 0.45, 0.40])
# Per-sample ILC-rich TLS proportion per niche
dist_mat = []
sample_names = []
for sid in sorted(niche_dist_samples.keys()):
    counts = niche_dist_samples[sid]
    total = sum(counts) + 1e-8
    dist_mat.append([c/total for c in counts])
    sample_names.append(sid)
dist_mat = np.array(dist_mat)
im_b = ax_b.imshow(dist_mat.T, aspect="auto", cmap="YlOrRd")
ax_b.set_xticks(range(len(sample_names)))
ax_b.set_xticklabels([s[:12] for s in sample_names], fontsize=3.5, rotation=90)
ax_b.set_yticks(range(K)); ax_b.set_yticklabels([f"N{n}" for n in range(K)], fontsize=5.5)
ax_b.set_title(f"b  ILC-rich TLS proportion per niche ({len(sample_names)} samples)", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_b, ax=ax_b, shrink=0.7)

# Panel C: Receptor detection heatmap (bottom-left)
df_rcp = df_all[df_all["type"]=="receptor"]
df_show = df_rcp[df_rcp["candidate"]].head(25)
if len(df_show) < 5: df_show = df_rcp.nlargest(25, "median_delta")
df_show = df_show.sort_values("median_delta", ascending=False)
rcp_mat = df_show[["pct_other","pct_ILC"]].values
ax_c = fig.add_axes([0.05, 0.06, 0.40, 0.40])
im_c = ax_c.imshow(rcp_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.6)
ax_c.set_xticks([0, 1]); ax_c.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_c.set_yticks(range(len(df_show))); ax_c.set_yticklabels(df_show["gene"].values, fontsize=5.5, fontstyle="italic")
ax_c.set_title(f"c  Receptor detection rate ({len(df_show)} candidates)", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_c, ax=ax_c, shrink=0.7)

# Panel D: Ligand detection heatmap (bottom-right)
df_lig = df_all[df_all["type"]=="ligand"]
df_lig_show = df_lig[df_lig["candidate"]].head(15)
if len(df_lig_show) < 3: df_lig_show = df_lig.nlargest(15, "median_delta")
df_lig_show = df_lig_show.sort_values("median_delta", ascending=False)
lig_mat = df_lig_show[["pct_other","pct_ILC"]].values
ax_d = fig.add_axes([0.52, 0.06, 0.20, 0.40])
im_d = ax_d.imshow(lig_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.6)
ax_d.set_xticks([0, 1]); ax_d.set_xticklabels(["other","ILC-rich"], fontsize=6)
ax_d.set_yticks(range(len(df_lig_show))); ax_d.set_yticklabels(df_lig_show["gene"].values, fontsize=5.5, fontstyle="italic")
ax_d.set_title(f"d  Ligand detection", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_d, ax=ax_d, shrink=0.7)

# Panel E: Category legend for receptors
ax_e = fig.add_axes([0.78, 0.06, 0.18, 0.40])
ax_e.axis("off")
cats_in_rcp = set(df_show["category"].values)
for j, (cat, color) in enumerate(CAT_COLORS.items()):
    if cat in cats_in_rcp:
        n_cat = sum(1 for _, r in df_show.iterrows() if r["category"]==cat)
        ax_e.text(0, 1.0 - j*0.08, f"{cat.split('(')[0].strip()} ({n_cat})", fontsize=5.5, color=color, fontweight="bold", transform=ax_e.transAxes)

fig.savefig(OUT / "fig_heatmap_final.svg", bbox_inches="tight")
fig.savefig(OUT / "fig_heatmap_final.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_heatmap_final.tiff", dpi=600, bbox_inches="tight")
plt.close()
print(f"Figure: {OUT}/fig_heatmap_final.{{svg,pdf,tiff}}")
print(f"Candidate receptors: {df_rcp['candidate'].sum()}, ligands: {df_lig['candidate'].sum()}")
