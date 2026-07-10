"""Receptor enrichment — all GBM samples with ILC-dominant TLS (Nature lollipop chart)"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

# Config
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)", "Inhib (GABA/Gly)", "Cholinergic (ACh)", "DA/NE", "Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna():
        gene_cat[g] = cat_order[idx]
cat_colors = {"Excit (Glutamate)":"#C24B3C","Inhib (GABA/Gly)":"#358554",
              "Cholinergic (ACh)":"#3575A3","DA/NE":"#764E9F","Serotonin (5-HT)":"#BF5A2E"}
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"),
     Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"),
     Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# === Per-sample per-gene FC ===
gene_fcs = {}

for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue

        tls = pd.read_csv(tls_csv)
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05, "values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        tls_mask = (tls["TLS.region"] == "TLS").values
        dominant = np.array(ct)[q05.argmax(axis=1)]

        ilc_dom = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES:
            idx = ct.index(c)
            p75 = np.percentile(q05[:, idx], 75)
            ilc_dom |= tls_mask & (dominant == c) & (q05[:, idx] >= p75)
        if ilc_dom.sum() < 3: continue

        ge = adata[:, adata.var["feature_types"] == "Gene Expression"] if "feature_types" in adata.var else adata
        expr = ge.X.toarray() if hasattr(ge.X, "toarray") else ge.X
        vn = ge.var_names.values

        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)
            im = np.log1p(expr[ilc_dom, gidx]).mean()
            cm = np.log1p(expr[~tls_mask, gidx]).mean()
            fc = im / cm if cm > 0 else np.nan
            gene_fcs.setdefault(g, []).append(fc)

# === Aggregate ===
rows = []
for g, fcs in gene_fcs.items():
    arr = np.array([x for x in fcs if not np.isnan(x)])
    if len(arr) < 5: continue
    median_fc = np.median(arr)
    log2fc = np.log2(median_fc) if median_fc > 0 else -np.inf
    try:
        _, pval = wilcoxon(arr - 1.0, alternative="two-sided")
    except:
        pval = 1.0
    rows.append({"gene": g, "n_samples": len(arr), "median_FC": median_fc,
                 "log2FC": log2fc, "pvalue": pval, "category": gene_cat.get(g, "?")})

df = pd.DataFrame(rows)
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("log2FC")
df.to_csv(Path(r"E:/GBM/results") / "receptor_gbm_all_per_gene.csv", index=False)

df_exp = df[df["log2FC"] > -10]
df_zero = df[df["log2FC"] <= -10]
n_max = df["n_samples"].max()

print(f"GBM samples with >=3 ILC-dominant TLS spots: {n_max}")
print(f"Genes: {len(df)} total, {len(df_exp)} expressed, {len(df_zero)} zero-expression")
print(f"FDR<0.1: {(df['fdr']<0.1).sum()} (enriched: {(df_exp['fdr']<0.1).sum()}, depleted: {(df_zero['fdr']<0.1).sum()})")

print("\n=== ENRICHED ===")
for _, r in df_exp[df_exp["fdr"] < 0.1].iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} log2FC={r['log2FC']:+6.3f}  FC={r['median_FC']:.2f}  FDR={r['fdr']:.4f}  n={r['n_samples']}  {s}")

print("\n=== 5-HT (zero expression) ===")
for _, r in df_zero[df_zero["category"].str.contains("5-HT")].iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*" if r["fdr"]<0.1 else ""
    print(f"  {r['gene']:8s} FC=0  FDR={r['fdr']:.4f}  n={r['n_samples']}  {s}")

# === FIGURE ===
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,
    "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.5,"legend.frameon":False})

fig, ax = plt.subplots(figsize=(5.5, 6.5))
N = len(df)

for i, (_, r) in enumerate(df.iterrows()):
    cat = r["category"]; color = cat_colors.get(cat, "gray")
    if r["log2FC"] <= -10:
        ax.axvline(x=-1.8, ymin=(i-0.35)/N, ymax=(i+0.35)/N, color="#cccccc", linewidth=0.5)
        ax.text(-1.9, i, r["gene"], ha="right", va="center", fontsize=5.2, fontstyle="italic", color="#999999")
        if r["fdr"] < 0.1:
            ax.scatter(-1.6, i, s=15, c="#cccccc", edgecolors="#999999", linewidths=0.4, zorder=5)
    else:
        xi = r["log2FC"]; alpha = 0.15 if r["fdr"]>=0.1 else 0.6
        ax.plot([0,xi],[i,i],color=color,linewidth=1.2 if alpha>0.3 else 0.5,alpha=alpha,zorder=3,solid_capstyle="round")
        is_sig = r["fdr"]<0.1
        ax.scatter(xi,i,s=32 if is_sig else 16,c=color,alpha=0.9 if is_sig else 0.3,
                   edgecolors="black" if is_sig else "none",linewidths=0.4,zorder=5)
        ax.text(xi-0.05,i,r["gene"],ha="right",va="center",fontsize=6,fontstyle="italic")

ax.axvline(x=0,color="black",linewidth=0.6,zorder=2)
ax.set_yticks([])
ax.set_xlabel("log2(ILC-TLS / non-TLS)",fontsize=8)
ax.set_title(f"Neurotransmitter receptors: ILC-dominant TLS vs non-TLS\nGBM only, n={n_max} samples",fontsize=9,fontweight="bold",loc="left")

ax.legend([Line2D([0],[0],marker="o",color="w",markerfacecolor="gray",markersize=6,markeredgecolor="black",markeredgewidth=0.4),
           Line2D([0],[0],marker="o",color="w",markerfacecolor="gray",markersize=3.5),
           Line2D([0],[0],marker="o",color="w",markerfacecolor="#cccccc",markersize=4,markeredgecolor="#999999",markeredgewidth=0.4)],
          ["FDR < 0.1","NS","Not expressed"],fontsize=5.5,loc="lower right",handletextpad=0.5)

for j, cat in enumerate(cat_order):
    n_ex = sum(1 for _,r in df_exp.iterrows() if r["category"]==cat)
    if n_ex > 0:
        ax.text(1.02, len(df_exp)+2-j*2.5, f"{cat} ({n_ex})", fontsize=5.5, color=cat_colors[cat],
                va="top",fontweight="bold",transform=ax.get_yaxis_transform())

fig.tight_layout(pad=0.5)
OUT = Path(r"E:\GBM\results")
for fmt in ["svg","pdf","tiff"]:
    fig.savefig(OUT / f"fig_receptor_gbm_final.{fmt}", bbox_inches="tight", dpi=600)
plt.close()
print(f"\nFigure: {OUT}/fig_receptor_gbm_final.{{svg,pdf,tiff}}")
