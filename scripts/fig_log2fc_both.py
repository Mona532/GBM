"""log2FC: ILC-enriched TLS vs non-TLS + niche-stratified comparison"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

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

def aggregate(gene_data):
    rows = []
    for g, vals in gene_data.items():
        arr = np.array(vals);
        if len(arr) < 5: continue
        _, p = wilcoxon(arr, alternative="two-sided")
        rows.append({"gene":g,"n_samples":len(arr),"median_log2FC":np.median(arr),"pvalue":p,"category":gene_cat.get(g,"?")})
    df = pd.DataFrame(rows)
    df["fdr"] = false_discovery_control(df["pvalue"].values)
    return df.sort_values("median_log2FC")

# ====== C1: ILC-enriched TLS vs non-TLS ======
g1 = {}
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
        if tls_mask.sum() < 5: continue
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}
        ilc_high = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx[c]] >= GLOBAL_THRESH[c])
        ilc_enr = ilc_high & tls_mask
        if ilc_enr.sum() < 3: continue
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_idx = [list(vn_all).index(g) for g in rx_present]
        rx_dense = ge[:, rx_idx].X.toarray() if hasattr(ge[:, rx_idx].X,"toarray") else np.asarray(ge[:, rx_idx].X)
        expr_norm = np.log1p(rx_dense / (lib_size[:,None]/10000))
        ie_m = expr_norm[ilc_enr].mean(axis=0)
        ne_m = expr_norm[~tls_mask].mean(axis=0)
        for gi, g in enumerate(rx_present):
            if ie_m[gi] > 0 and ne_m[gi] > 0:
                g1.setdefault(g, []).append(np.log2(ie_m[gi]/ne_m[gi]))
df1 = aggregate(g1)
df1.to_csv(OUT / "receptor_log2fc_ilc_vs_nontls.csv", index=False)
print(f"C1 (ILC vs non-TLS): {len(df1)} total, {(df1['fdr']<0.1).sum()} FDR<0.1")

# ====== C2: Niche-stratified ILC-enriched vs other TLS ======
g2 = {}
nd = 0
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
        ilc_idx = {c: ct_names.index(c) for c in ILC_TYPES}
        coords = adata.obsm["spatial"]
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_idx = [list(vn_all).index(g) for g in rx_present]
        rx_dense = ge[:, rx_idx].X.toarray() if hasattr(ge[:, rx_idx].X,"toarray") else np.asarray(ge[:, rx_idx].X)
        expr_norm = np.log1p(rx_dense / (lib_size[:,None]/10000))

        obs_df = pd.DataFrame(index=adata.obs_names); obs_df["sample"] = sd.name
        a = ad.AnnData(X=q05.astype(np.float32), obs=obs_df, obsm={"spatial": coords})
        a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
        try:
            sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
            cc.gr.aggregate_neighbors(a, n_layers=1, use_rep=None, out_key="X_cc", sample_key="sample")
            autok = cc.tl.ClusterAutoK(n_clusters=(5,15), max_runs=3)
            autok.fit(a, use_rep="X_cc")
            niches = autok.predict(a, use_rep="X_cc")
        except: continue

        for niche_id in sorted(set(niches)):
            niche_mask = niches == niche_id
            tls_in_niche = tls_mask & niche_mask
            if tls_in_niche.sum() < 10: continue
            tls_spot_idx = np.where(tls_in_niche)[0]
            ilc_enr = np.zeros(len(tls_spot_idx), dtype=bool)
            for c in ILC_TYPES: ilc_enr |= (q05[tls_spot_idx, ilc_idx[c]] >= GLOBAL_THRESH[c])
            other_tls = ~ilc_enr
            if ilc_enr.sum() < 3 or other_tls.sum() < 3: continue
            ie_m = expr_norm[tls_spot_idx[ilc_enr]].mean(axis=0)
            oe_m = expr_norm[tls_spot_idx[other_tls]].mean(axis=0)
            for gi, g in enumerate(rx_present):
                if ie_m[gi] > 0 and oe_m[gi] > 0:
                    g2.setdefault(g, []).append(np.log2(ie_m[gi]/oe_m[gi]))
        nd += 1
        if nd % 20 == 0: print(f"  niche: {nd}")
print(f"  niche: {nd} samples total")
df2 = aggregate(g2)
df2.to_csv(OUT / "receptor_log2fc_niche_stratified.csv", index=False)
print(f"C2 (niche-stratified): {len(df2)} total, {(df2['fdr']<0.1).sum()} FDR<0.1")

# ====== Figure ======
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))

for ax, df_in, title, label in [
    (ax1, df1, "ILC-enriched TLS vs non-TLS", "a"),
    (ax2, df2, "Niche-stratified\n(ILC-enriched vs other TLS)", "b")]:
    df_pos = df_in[(df_in["median_log2FC"] > 0.2) & (df_in["fdr"] < 0.1)]
    if len(df_pos) == 0:
        df_pos = df_in.nlargest(15, "median_log2FC")
    df_plot = df_pos.sort_values("median_log2FC", ascending=True)
    N = len(df_plot)
    y = np.arange(N)
    for i, (_, r) in enumerate(df_plot.iterrows()):
        c = CAT.get(r["category"], "#888")
        is_sig = r["fdr"] < 0.1
        ax.plot([0, r["median_log2FC"]], [i, i], color=c, linewidth=1.5, alpha=0.6 if is_sig else 0.2, solid_capstyle="round", zorder=3)
        ax.scatter(r["median_log2FC"], i, s=45 if is_sig else 25, c=c, alpha=0.9 if is_sig else 0.3,
                   edgecolors="black" if is_sig else "none", linewidths=0.4, zorder=5)
    ax.set_yticks(range(N))
    ax.set_yticklabels(df_plot["gene"].values, fontsize=6.5, fontstyle="italic")
    ax.axvline(x=0, color="black", linewidth=0.4)
    ax.set_xlabel("median log2FC", fontsize=7)
    ax.set_title(f"{label}  {title}", fontsize=8, fontweight="bold", loc="left")
    if title.startswith("Niche"):
        ax.text(0.98, 0.02, "All FDR > 0.1", transform=ax.transAxes, fontsize=6, color="#c44e52", ha="right", fontweight="bold")

cats_in = set(df1[df1["fdr"]<0.1]["category"].values) | set(df2[df2["fdr"]<0.1]["category"].values)
fig.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in],
           fontsize=6, loc="lower center", ncol=5, borderpad=0.3)
fig.tight_layout(pad=0.8, rect=[0, 0.06, 1, 1])
save_pub(fig, OUT / "fig_receptor_log2fc_comparison")
plt.close()
print(f"Figure: {OUT}/fig_receptor_log2fc_comparison.{{svg,pdf,tiff}}")
