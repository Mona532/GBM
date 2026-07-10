"""
Final main figure: CellCharter niches as spatial context, sample-level receptor comparison.
ILC-enriched TLS vs other TLS, 6 positive criteria.
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import gc, warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})
def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

CAT = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): gene_cat[g] = list(CAT.keys())[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]
GLOBAL_THRESH = {"ILC1":1.034,"ILC2":1.0,"ILC3":1.035}
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results")
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# ====== Step 1: Global CellCharter niche (z-scored q05) ======
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
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_idx = [list(vn_all).index(g) for g in rx_present]
        rx_raw = ge[:, rx_idx].X.toarray() if hasattr(ge[:, rx_idx].X,"toarray") else np.asarray(ge[:, rx_idx].X)
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))
        sample_data.append((sd.name, tls_mask, q05, expr_norm, rx_raw, rx_present, ct_names, ilc_idx_map, coords))

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

# ====== Step 2: Per-sample (aggregated across niches) ILC-enriched vs other TLS ======
gene_sample = {}        # gene -> [{sample, delta, l2fc, pct_ie, pct_ot}]
niche_distribution = {} # sample -> niche distribution of ILC-enriched TLS

offset = 0
for sid, tls_mask, q05, expr_norm, rx_raw, rx_present, ct_names, ilc_idx_map, coords in sample_data:
    n_spots = len(tls_mask)
    niches = niche_all[offset:offset+n_spots]
    offset += n_spots

    # ILC-enriched TLS (across all niches)
    ilc_high = np.zeros(n_spots, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
    ilc_enr = ilc_high & tls_mask
    other_tls = tls_mask & ~ilc_enr
    if ilc_enr.sum() < 3 or other_tls.sum() < 3: continue

    # ILC-enriched TLS distribution across niches
    niche_distribution[sid] = {}
    for n in range(K):
        n_in_niche = (niches == n) & ilc_enr
        niche_distribution[sid][n] = n_in_niche.sum()

    # Per-gene comparison: ILC-enriched vs other TLS (sample-level, all niches merged)
    ie_expr = expr_norm[ilc_enr]; ot_expr = expr_norm[other_tls]
    for gi, g in enumerate(rx_present):
        ie_m = ie_expr[:, gi].mean(); ot_m = ot_expr[:, gi].mean()
        delta = ie_m - ot_m; l2fc = np.log2((ie_m+1e-6)/(ot_m+1e-6))
        pct_ie = (rx_raw[ilc_enr, gi] > 0).mean()
        pct_ot = (rx_raw[other_tls, gi] > 0).mean()
        gene_sample.setdefault(g, []).append({"sample":sid,"delta":delta,"log2fc":l2fc,"pct_ie":pct_ie,"pct_ot":pct_ot})

# ====== Step 3: Positive receptor criteria ======
rows = []
for g, svals in gene_sample.items():
    n = len(svals); deltas = np.array([v["delta"] for v in svals])
    pct_ies = np.array([v["pct_ie"] for v in svals]); pct_ots = np.array([v["pct_ot"] for v in svals])
    l2fcs = np.array([v["log2fc"] for v in svals])
    _, p_delta = wilcoxon(deltas, alternative="two-sided")
    rows.append({"gene":g,"n":n,"median_delta":np.median(deltas),"median_log2FC":np.median(l2fcs),
        "pct_ILC_TLS":np.median(pct_ies),"pct_other_TLS":np.median(pct_ots),
        "prop_positive":(deltas>0).mean(),"p_delta":p_delta,"category":gene_cat.get(g,"?")})

df_all = pd.DataFrame(rows)
df_all = df_all[df_all["p_delta"].notna() & (df_all["p_delta"]>0) & (df_all["p_delta"]<1)]
df_all["fdr_delta"] = false_discovery_control(df_all["p_delta"].values)

POSITIVE = ((df_all["n"]>=20) & (df_all["pct_ILC_TLS"]>=0.05) & (df_all["median_delta"]>0) &
            (df_all["fdr_delta"]<0.1) & (df_all["prop_positive"]>=0.6) & (df_all["pct_ILC_TLS"]>df_all["pct_other_TLS"]))
df_all["positive"] = POSITIVE
df_pos = df_all[POSITIVE].sort_values("median_delta", ascending=False)
df_all.to_csv(OUT / "receptor_sample_level_final.csv", index=False)
print(f"Positive receptors: {len(df_pos)}")
for _, r in df_pos.iterrows():
    print(f"  {r['gene']:8s} {r['category']:22s} delta={r['median_delta']:+.4f} pct_IE={r['pct_ILC_TLS']:.3f} pct_OT={r['pct_other_TLS']:.3f} prop+={r['prop_positive']:.2f} n={r['n']}")

# ====== Figure: 4-panel ======
fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(10, 8))

# --- Panel A: Niche spatial map (demo sample) ---
demo_sid, demo_mask, demo_q05, demo_expr, demo_raw, demo_rx, demo_ct, demo_ilc, demo_coords = sample_data[0]
off = sum(s[1].sum() for s in sample_data[:0])
demo_niches = niche_all[off:off+len(demo_mask)]
cmap = ListedColormap(plt.cm.tab10.colors[:K])
ax_a.scatter(demo_coords[:,0], demo_coords[:,1], c=demo_niches.astype(int), s=0.5, cmap=cmap, rasterized=True)
ax_a.set_title(f"a  CellCharter niches ({K} global)\n{demo_sid}", fontsize=7.5, fontweight="bold", loc="left")
ax_a.axis("off")

# --- Panel B: ILC-enriched TLS distribution across niches ---
niche_counts = {n: 0 for n in range(K)}
for sid, dist in niche_distribution.items():
    for n, c in dist.items(): niche_counts[n] += c
niche_df = pd.DataFrame({"niche": list(niche_counts.keys()), "count": list(niche_counts.values())})
niche_df = niche_df.sort_values("count", ascending=False)
bars = ax_b.bar(range(K), niche_df["count"].values, color=[plt.cm.tab10.colors[i % 10] for i in niche_df["niche"].values], alpha=0.6)
ax_b.set_xticks(range(K))
ax_b.set_xticklabels([f"N{n}" for n in niche_df["niche"].values], fontsize=6)
ax_b.set_ylabel("ILC-enriched TLS spots", fontsize=7)
ax_b.set_title(f"b  ILC-enriched TLS across niches ({len(niche_distribution)} samples)", fontsize=7.5, fontweight="bold", loc="left")

# --- Panel C: Receptor lollipop ---
df_plot = df_pos.head(20) if len(df_pos) > 0 else df_all.nlargest(20, "median_delta").query("median_delta > 0")
df_plot = df_plot.sort_values("median_delta", ascending=True)
y = np.arange(len(df_plot))
for i, (_, r) in enumerate(df_plot.iterrows()):
    c = CAT.get(r["category"], "#888")
    is_pos = r["positive"]
    sz = 15 + max(r["pct_ILC_TLS"], 0.01) * 150
    ax_c.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.6 if is_pos else 0.15, solid_capstyle="round", zorder=3)
    ax_c.scatter(r["median_delta"], i, s=sz, c=c, alpha=0.9 if is_pos else 0.3, edgecolors="black" if is_pos else "none", linewidths=0.4, zorder=5)
ax_c.set_yticks(range(len(df_plot)))
ax_c.set_yticklabels(df_plot["gene"].values, fontsize=6.5, fontstyle="italic")
ax_c.axvline(x=0, color="black", linewidth=0.4)
ax_c.set_xlabel("median delta (ILC-enriched - other TLS)", fontsize=7)
ax_c.set_title(f"c  Receptor expression ({len(df_pos)} positive, 6 criteria)", fontsize=7.5, fontweight="bold", loc="left")

# --- Panel D: Paired for top 4 ---
top_genes = df_pos.head(4) if len(df_pos) >= 4 else df_all.nlargest(4, "median_delta")
top4g = [g for g in top_genes["gene"].values if g in gene_sample and len(gene_sample[g]) > 5][:4]
for gi, g in enumerate(top4g):
    sv = gene_sample[g]
    ie_v = np.array([v["pct_ie"] for v in sv]); ot_v = np.array([v["pct_ot"] for v in sv])
    for i in range(min(len(ie_v), 100)):
        ax_d.plot([gi-0.25, gi+0.25], [ot_v[i], ie_v[i]], color="#ccc", linewidth=0.2, alpha=0.3, zorder=2)
    ax_d.scatter(np.full(len(ot_v), gi-0.25), ot_v, s=6, c="#4c72b0", alpha=0.3, zorder=3)
    ax_d.scatter(np.full(len(ie_v), gi+0.25), ie_v, s=6, c="#c44e52", alpha=0.3, zorder=3)
    ax_d.scatter(gi-0.25, np.median(ot_v), s=40, c="#4c72b0", edgecolors="black", linewidths=0.5, zorder=6)
    ax_d.scatter(gi+0.25, np.median(ie_v), s=40, c="#c44e52", edgecolors="black", linewidths=0.5, zorder=6)
ax_d.set_xticks(range(len(top4g)))
ax_d.set_xticklabels(top4g, fontsize=7, fontstyle="italic")
ax_d.set_ylabel("pct_expr", fontsize=7)
ax_d.set_title("d  Paired per-sample: other TLS vs ILC-enriched", fontsize=7.5, fontweight="bold", loc="left")
ax_d.legend([Line2D([0],[0],marker='o',c='w',markerfacecolor='#4c72b0',markersize=8),
             Line2D([0],[0],marker='o',c='w',markerfacecolor='#c44e52',markersize=8)],
            ["other TLS","ILC-enriched"], fontsize=6, loc="upper left")

fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_final_main_4panel")
plt.close()
print(f"Figure: {OUT}/fig_final_main_4panel.{{svg,pdf,tiff}}")
