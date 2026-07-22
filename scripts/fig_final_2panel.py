"""Final 2-panel figure: (a) ILC-enriched vs non-TLS log2FC, (b) global CellCharter niche-stratified delta"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})
def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

CAT = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}
OUT = Path(r"E:/GBM/results")

# ====== Panel A: ILC-enriched TLS vs non-TLS log2FC (from C1, correctly computed earlier) ======
# Re-run C1 cleanly
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = list(CAT.keys())
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): gene_cat[g] = cat_order[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]
GLOBAL_THRESH = {"ILC1":1.034,"ILC2":1.0,"ILC3":1.035}
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}

import anndata as ad
from scipy.stats import wilcoxon, false_discovery_control

DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

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

rows1 = []
for g, vals in g1.items():
    arr = np.array(vals)
    if len(arr) < 5: continue
    _, p = wilcoxon(arr, alternative="two-sided")
    rows1.append({"gene":g,"n_samples":len(arr),"median_log2FC":np.median(arr),"pvalue":p,"category":gene_cat.get(g,"?")})
df1 = pd.DataFrame(rows1)
df1["fdr"] = false_discovery_control(df1["pvalue"].values)
df1.to_csv(OUT / "fig_panel_a_log2fc.csv", index=False)
print(f"Panel A: {len(df1)} genes, {(df1['fdr']<0.1).sum()} FDR<0.1")

# ====== Panel B: Global CellCharter niche-stratified delta (from official results) ======
df2 = pd.read_csv(OUT / "receptor_cellcharter_final.csv")
print(f"Panel B: {len(df2)} genes, {(df2['fdr']<0.1).sum()} FDR<0.1")

# ====== Figure ======
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))

# Panel A: log2FC, only enriched FDR<0.1
df_a = df1[(df1["median_log2FC"] > 1.5) & (df1["fdr"] < 0.1)].sort_values("median_log2FC", ascending=True)
N = len(df_a)
y = np.arange(N)
for i, (_, r) in enumerate(df_a.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax1.plot([0, r["median_log2FC"]], [i, i], color=c, linewidth=1.5, alpha=0.6, solid_capstyle="round", zorder=3)
    ax1.scatter(r["median_log2FC"], i, s=50, c=c, alpha=0.9, edgecolors="black", linewidths=0.4, zorder=5)
ax1.set_yticks(range(N))
ax1.set_yticklabels(df_a["gene"].values, fontsize=7, fontstyle="italic")
ax1.axvline(x=0, color="black", linewidth=0.4)
ax1.set_xlabel("median log2FC", fontsize=7)
ax1.set_title(f"a  ILC-enriched TLS vs non-TLS\n(n={df_a['n_samples'].max()}, {N} genes FDR<0.1)", fontsize=8, fontweight="bold", loc="left")

# Panel B: niche-stratified delta, all genes (null result)
df_b = df2.nlargest(20, "median_delta").sort_values("median_delta", ascending=True)
N2 = len(df_b)
y2 = np.arange(N2)
for i, (_, r) in enumerate(df_b.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax2.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.2, solid_capstyle="round", zorder=3)
    ax2.scatter(r["median_delta"], i, s=40, c=c, alpha=0.5, edgecolors="#999", linewidths=0.3, zorder=5)
ax2.set_yticks(range(N2))
ax2.set_yticklabels(df_b["gene"].values, fontsize=7, fontstyle="italic")
ax2.axvline(x=0, color="black", linewidth=0.4)
ax2.set_xlabel("median delta (log-expr)", fontsize=7)
ax2.set_title(f"b  CellCharter niche-stratified, global 13 niches\n(n={df_b['n_samples'].max()}, 0 genes FDR<0.1)", fontsize=8, fontweight="bold", loc="left")
ax2.text(0.98, 0.02, "All FDR > 0.1", transform=ax2.transAxes, fontsize=7, color="#c44e52", ha="right", fontweight="bold")

# Shared category legend
cats_in = set(df_a["category"].values) | set(df_b["category"].values)
fig.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in],
           fontsize=6, loc="lower center", ncol=5, borderpad=0.3)
fig.tight_layout(pad=0.8, rect=[0, 0.06, 1, 1])
save_pub(fig, OUT / "fig_final_2panel")
plt.close()
print(f"Figure: {OUT}/fig_final_2panel.{{svg,pdf,tiff}}")
