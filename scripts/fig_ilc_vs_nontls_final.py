"""ILC-enriched TLS vs non-TLS — 6 positive receptor criteria, Nature figure"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
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

gene_sample = {}
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
        rx_raw = ge[:, rx_idx].X.toarray() if hasattr(ge[:, rx_idx].X,"toarray") else np.asarray(ge[:, rx_idx].X)
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))
        ie = expr_norm[ilc_enr]; ne = expr_norm[~tls_mask]
        for gi, g in enumerate(rx_present):
            ie_m = ie[:, gi].mean(); ne_m = ne[:, gi].mean()
            delta = ie_m - ne_m; l2fc = np.log2((ie_m+1e-6)/(ne_m+1e-6))
            pct_ie = (rx_raw[ilc_enr, gi] > 0).mean()
            pct_ne = (rx_raw[~tls_mask, gi] > 0).mean()
            gene_sample.setdefault(g, []).append({"sample":sd.name,"delta":delta,"log2fc":l2fc,"pct_ie":pct_ie,"pct_ne":pct_ne})

rows = []
for g, svals in gene_sample.items():
    n = len(svals)
    deltas = np.array([v["delta"] for v in svals])
    pct_ies = np.array([v["pct_ie"] for v in svals]); pct_nes = np.array([v["pct_ne"] for v in svals])
    l2fcs = np.array([v["log2fc"] for v in svals])
    _, p_delta = wilcoxon(deltas, alternative="two-sided")
    rows.append({"gene":g,"n":n,"median_delta":np.median(deltas),"median_log2FC":np.median(l2fcs),
        "pct_ILC_TLS":np.median(pct_ies),"pct_nonTLS":np.median(pct_nes),
        "prop_positive":(deltas>0).mean(),"p_delta":p_delta,"category":gene_cat.get(g,"?")})

df_all = pd.DataFrame(rows)
df_all = df_all[df_all["p_delta"].notna() & (df_all["p_delta"]>0) & (df_all["p_delta"]<1)]
df_all["fdr_delta"] = false_discovery_control(df_all["p_delta"].values)

POSITIVE = ((df_all["n"]>=20) & (df_all["pct_ILC_TLS"]>=0.05) & (df_all["median_delta"]>0) &
            (df_all["fdr_delta"]<0.1) & (df_all["prop_positive"]>=0.6) & (df_all["pct_ILC_TLS"]>df_all["pct_nonTLS"]))
df_all["positive"] = POSITIVE
df_pos = df_all[POSITIVE].sort_values("median_delta", ascending=False)
df_all.to_csv(OUT / "receptor_ilc_vs_nontls_positive.csv", index=False)
print(f"Positive receptors: {len(df_pos)}")
for _, r in df_pos.iterrows():
    print(f"  {r['gene']:8s} {r['category']:22s} delta={r['median_delta']:+.4f} pct_IE={r['pct_ILC_TLS']:.3f} pct_NE={r['pct_nonTLS']:.3f} prop+={r['prop_positive']:.2f} n={r['n']} FDR={r['fdr_delta']:.4f}")

df_plot = df_pos.head(20) if len(df_pos) > 0 else df_all.nlargest(20, "median_delta")
df_plot = df_plot.sort_values("median_delta", ascending=True)
fig, ax = plt.subplots(figsize=(5, len(df_plot)*0.33+1))
y = np.arange(len(df_plot))
for i, (_, r) in enumerate(df_plot.iterrows()):
    c = CAT.get(r["category"], "#888")
    is_pos = r["positive"]
    sz = 15 + max(r["pct_ILC_TLS"], 0.01) * 150
    ax.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.6 if is_pos else 0.15, solid_capstyle="round", zorder=3)
    ax.scatter(r["median_delta"], i, s=sz, c=c, alpha=0.9 if is_pos else 0.3, edgecolors="black" if is_pos else "none", linewidths=0.4, zorder=5)
    if is_pos: ax.text(r["median_delta"] + 0.002, i, f"n={r['n']}", fontsize=5.5, color="#888", va="center")
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot["gene"].values, fontsize=7, fontstyle="italic")
ax.axvline(x=0, color="black", linewidth=0.4)
ax.set_xlabel("median delta (ILC-enriched TLS - non-TLS)", fontsize=8)
ax.set_title(f"ILC-enriched TLS vs non-TLS ({len(df_pos)} positive, 6 criteria)", fontsize=9, fontweight="bold", loc="left")
cats_in = set(df_plot["category"].values)
ax.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in], fontsize=6, loc="lower right", borderpad=0.3)
fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_ilc_vs_nontls_positive")
plt.close()
print(f"Figure: {OUT}/fig_ilc_vs_nontls_positive.{{svg,pdf,tiff}}")
