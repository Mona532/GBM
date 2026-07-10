"""
Final complete analysis: ILC-rich TLS microenvironment — receptor & ligand spatial signals
CellCharter for spatial context, sample-level receptor/ligand comparison
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import lil_matrix
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
RECEPTOR_GENES = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): RECEPTOR_GENES[g] = list(CAT.keys())[idx]

# Ligand-receptor pairs
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

# ====== Step 1: Global CellCharter ======
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
        sample_data.append((sd.name, tls_mask, q05, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords))

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

# ====== Step 2: ILC-rich TLS definition + receptor/ligand comparison ======
gene_sample = {}     # gene -> [{sample, delta, l2fc, pct_ilc, pct_other}]
ligand_sample = {}   # ligand -> [{sample, delta, pct_ilc_nhood, pct_other_nhood}]
niche_dist = {}
offset = 0

for sid, tls_mask, q05, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords in sample_data:
    n = len(tls_mask); niches = niche_all[offset:offset+n]; offset += n

    # ILC-rich TLS = any ILC >= max(global P75, 1.0)
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
    ilc_rich = ilc_high & tls_mask
    other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue

    niche_dist[sid] = {ni: int((niches == ni).sum()) for ni in range(K)}

    # Gene expression comparison
    ie = expr_norm[ilc_rich]; oe = expr_norm[other_tls]
    for gi, g in enumerate(genes_present):
        delta = ie[:, gi].mean() - oe[:, gi].mean()
        l2fc = np.log2((ie[:, gi].mean() + 1e-6) / (oe[:, gi].mean() + 1e-6))
        pct_ilc = (rx_raw[ilc_rich, gi] > 0).mean()
        pct_other = (rx_raw[other_tls, gi] > 0).mean()
        entry = {"sample":sid,"delta":delta,"log2fc":l2fc,"pct_ilc":pct_ilc,"pct_other":pct_other}
        if g in RECEPTOR_GENES:
            gene_sample.setdefault(g, []).append(entry)
        if g in LIGANDS:
            ligand_sample.setdefault(g, []).append(entry)

# ====== Step 3: Positive receptor/ligand criteria ======
def compute_positive(data_dict, gene_type="receptor"):
    rows = []
    for g, svals in data_dict.items():
        n = len(svals); deltas = np.array([v["delta"] for v in svals])
        pct_ilc = np.array([v["pct_ilc"] for v in svals]); pct_other = np.array([v["pct_other"] for v in svals])
        _, p_delta = wilcoxon(deltas, alternative="two-sided")
        rows.append({"gene":g,"n":n,"median_delta":np.median(deltas),
            "pct_ILC":np.median(pct_ilc),"pct_other":np.median(pct_other),
            "prop_positive":(deltas>0).mean(),"p_delta":p_delta,
            "category":RECEPTOR_GENES.get(g,"Ligand"),"type":gene_type})
    df = pd.DataFrame(rows)
    df = df[df["p_delta"].notna() & (df["p_delta"]>0) & (df["p_delta"]<1)]
    df["fdr"] = false_discovery_control(df["p_delta"].values)
    df["positive"] = ((df["n"]>=20)&(df["pct_ILC"]>=0.05)&(df["median_delta"]>0)&
                       (df["fdr"]<0.1)&(df["prop_positive"]>=0.6)&(df["pct_ILC"]>df["pct_other"]))
    df["candidate"] = ((df["n"]>=15)&(df["pct_ILC"]>=0.03)&(df["median_delta"]>0)&(df["prop_positive"]>=0.55))
    return df.sort_values("median_delta", ascending=False)

df_rcp = compute_positive(gene_sample, "receptor")
df_lig = compute_positive(ligand_sample, "ligand")
df_all = pd.concat([df_rcp, df_lig])
df_all.to_csv(OUT / "receptor_ligand_combined.csv", index=False)

n_pos = df_rcp["positive"].sum(); n_cand = df_rcp["candidate"].sum()
n_lig_pos = df_lig["positive"].sum(); n_lig_cand = df_lig["candidate"].sum()
print(f"Receptors: {n_pos} positive (strict), {n_cand} candidate")
print(f"Ligands: {n_lig_pos} positive (strict), {n_lig_cand} candidate")

if n_pos > 0:
    print("  Positive receptors:")
    for _, r in df_rcp[df_rcp["positive"]].iterrows():
        print(f"    {r['gene']:8s} delta={r['median_delta']:+.4f} pct_ILC={r['pct_ILC']:.3f} pct_OT={r['pct_other']:.3f} prop+={r['prop_positive']:.2f} n={r['n']}")
if n_cand > n_pos:
    print("  Additional candidate receptors:")
    cand_only = df_rcp[df_rcp["candidate"] & ~df_rcp["positive"]]
    for _, r in cand_only.head(10).iterrows():
        print(f"    {r['gene']:8s} delta={r['median_delta']:+.4f} pct_ILC={r['pct_ILC']:.3f} prop+={r['prop_positive']:.2f} n={r['n']}")

# ====== Figure ======
fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(11, 8.5))

# Panel A: Niche map
demo_sid, demo_mask, demo_q05, demo_expr, demo_raw, demo_genes, demo_ct, demo_ilc, demo_coords = sample_data[0]
off0 = sum(s[1].sum() for s in sample_data[:0])
demo_niches = niche_all[off0:off0+len(demo_mask)].astype(int)
cmap = ListedColormap(plt.cm.tab10.colors[:K])
ax_a.scatter(demo_coords[:,0], demo_coords[:,1], c=demo_niches, s=0.5, cmap=cmap, rasterized=True)
ax_a.set_title(f"a  CellCharter spatial niches (K={K})\n{demo_sid}", fontsize=8, fontweight="bold", loc="left")
ax_a.axis("off")

# Panel B: ILC-rich TLS niche distribution
niche_counts = {ni: 0 for ni in range(K)}
for sid, nd in niche_dist.items():
    for ni, c in nd.items(): niche_counts[ni] += c
niche_df2 = pd.DataFrame({"niche": list(niche_counts.keys()), "count": list(niche_counts.values())})
niche_df2 = niche_df2.sort_values("niche")
bars2 = ax_b.bar(range(K), niche_df2["count"].values, color=[plt.cm.tab10(i) for i in range(K)], alpha=0.6, edgecolor="black", linewidth=0.3)
ax_b.set_xticks(range(K))
ax_b.set_xticklabels([f"N{n}" for n in niche_df2["niche"].values], fontsize=6)
ax_b.set_ylabel("ILC-rich TLS spots", fontsize=7)
ax_b.set_title(f"b  ILC-rich TLS across niches ({len(niche_dist)} samples)", fontsize=8, fontweight="bold", loc="left")

# Panel C: Receptor + ligand lollipop (candidate level)
# Show candidates (less strict) since strict positive may be 0
df_show = df_all[df_all["candidate"]].head(25)
if len(df_show) == 0: df_show = df_all.nlargest(25, "median_delta")
df_show = df_show.sort_values("median_delta", ascending=True)
y = np.arange(len(df_show))
for i, (_, r) in enumerate(df_show.iterrows()):
    is_rcp = r["type"] == "receptor"
    c = CAT.get(r["category"], "#888") if is_rcp else "#888"
    marker = "o" if is_rcp else "s"
    sz = 15 + max(r["pct_ILC"], 0.01) * 150
    is_pos = r["positive"]
    ax_c.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.6 if is_pos else 0.2, solid_capstyle="round", zorder=3)
    ax_c.scatter(r["median_delta"], i, s=sz, c=c, marker=marker, alpha=0.9 if is_pos else 0.3,
                 edgecolors="black" if is_pos else "none", linewidths=0.4, zorder=5)
    label = f"{r['gene']}{'*' if is_rcp else '†'}"
    ax_c.text(r["median_delta"] + 0.002, i, label, fontsize=5.5, fontstyle="italic", va="center")
ax_c.set_yticks(range(len(df_show))); ax_c.set_yticklabels([])
ax_c.axvline(x=0, color="black", linewidth=0.4)
ax_c.set_xlabel("median delta (ILC-rich - other TLS)", fontsize=7)
ax_c.set_title(f"c  Receptor (*) & ligand (†) signals\n{len(df_show)} candidates ({n_pos} strict positive, {n_lig_pos} ligand)", fontsize=8, fontweight="bold", loc="left")

# Panel D: Top 4 paired
top4 = [g for _, r in df_show.head(4).iterrows() if (g := r["gene"]) in gene_sample or g in ligand_sample]
top4_data = {g: gene_sample.get(g, ligand_sample.get(g, [])) for g in top4}
top4g = [g for g in top4 if len(top4_data[g]) > 5][:4]
for gi, g in enumerate(top4g):
    sv = top4_data[g]; ie_v = np.array([v["pct_ilc"] for v in sv]); ot_v = np.array([v["pct_other"] for v in sv])
    for i in range(min(len(ie_v), 100)):
        ax_d.plot([gi-0.25, gi+0.25], [ot_v[i], ie_v[i]], color="#ccc", linewidth=0.2, alpha=0.3, zorder=2)
    ax_d.scatter(np.full(len(ot_v), gi-0.25), ot_v, s=6, c="#4c72b0", alpha=0.3, zorder=3)
    ax_d.scatter(np.full(len(ie_v), gi+0.25), ie_v, s=6, c="#c44e52", alpha=0.3, zorder=3)
    ax_d.scatter(gi-0.25, np.median(ot_v), s=45, c="#4c72b0", edgecolors="black", linewidths=0.5, zorder=6)
    ax_d.scatter(gi+0.25, np.median(ie_v), s=45, c="#c44e52", edgecolors="black", linewidths=0.5, zorder=6)
ax_d.set_xticks(range(len(top4g))); ax_d.set_xticklabels(top4g, fontsize=7, fontstyle="italic")
ax_d.set_ylabel("pct_expr", fontsize=7)
ax_d.set_title("d  Paired per-sample: other TLS vs ILC-rich TLS", fontsize=8, fontweight="bold", loc="left")
ax_d.legend([Line2D([0],[0],marker='o',c='w',markerfacecolor='#4c72b0',markersize=8),
             Line2D([0],[0],marker='o',c='w',markerfacecolor='#c44e52',markersize=8)],
            ["other TLS","ILC-rich TLS"], fontsize=6, loc="upper left")

fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_final_complete")
plt.close()
print(f"Figure: {OUT}/fig_final_complete.{{svg,pdf,tiff}}")
