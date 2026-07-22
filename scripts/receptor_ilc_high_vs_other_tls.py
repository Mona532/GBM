"""Receptor enrichment: ILC-high TLS vs other TLS (within-TLS comparison)"""
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

# ── Config ──
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)", "Inhib (GABA/Gly)", "Cholinergic (ACh)", "DA/NE", "Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna():
        gene_cat[g] = cat_order[idx]
cat_colors = {"Excit (Glutamate)":"#c44e52","Inhib (GABA/Gly)":"#55a868",
              "Cholinergic (ACh)":"#4c72b0","DA/NE":"#937860","Serotonin (5-HT)":"#ccb974"}
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"),
     Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"),
     Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# ── Step 1: Global TLS P75 per ILC ──
all_tls_ilc = {c: [] for c in ILC_TYPES}
for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() == 0: continue
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        for c in ILC_TYPES:
            all_tls_ilc[c].extend(q05[tls_mask, ct.index(c)].tolist())

GLOBAL_P75 = {c: np.percentile(all_tls_ilc[c], 75) for c in ILC_TYPES}
THRESH = {c: max(GLOBAL_P75[c], 1.0) for c in ILC_TYPES}
print(f"Thresholds: { {c: round(v,3) for c,v in THRESH.items()} }")

# ── Step 2: Per-sample per-gene FC: ILC-high TLS vs other TLS ──
gene_fcs = {}
n_ilc_samples = 0

for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 6: continue  # need both groups large enough
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}

        # ILC-high TLS mask
        ilc_mat = q05[:, [ilc_idx[c] for c in ILC_TYPES]]
        ilc_argmax = np.array(ILC_TYPES)[ilc_mat.argmax(axis=1)]
        ilc_high = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES:
            ilc_high |= tls_mask & (ilc_argmax == c) & (q05[:, ilc_idx[c]] >= THRESH[c])
        other_tls = tls_mask & ~ilc_high

        if ilc_high.sum() < 3 or other_tls.sum() < 3:
            continue
        n_ilc_samples += 1

        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values

        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)
            im = np.log1p(expr[ilc_high, gidx]).mean()
            cm = np.log1p(expr[other_tls, gidx]).mean()
            fc = im / cm if cm > 0 else np.nan
            gene_fcs.setdefault(g, []).append(fc)

# ── Aggregate ──
rows = []
for g, fcs in gene_fcs.items():
    arr = np.array([x for x in fcs if not np.isnan(x)])
    if len(arr) < 5: continue
    median_fc = np.median(arr)
    log2fc = np.log2(median_fc) if median_fc > 0 else -np.inf
    _, pval = wilcoxon(arr - 1.0, alternative="two-sided")
    rows.append({"gene": g, "n_samples": len(arr), "median_FC": median_fc,
                 "log2FC": log2fc, "pvalue": pval, "category": gene_cat.get(g, "?")})

df = pd.DataFrame(rows)
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("log2FC")
df.to_csv(Path(r"E:/GBM/results") / "receptor_ilc_high_vs_other_tls.csv", index=False)

df_exp = df[df["log2FC"] > -10]
n_max = df["n_samples"].max()
print(f"\nSamples with ILC-high TLS: {n_ilc_samples} (effective: {n_max})")
print(f"Genes: {len(df)}, expressed: {len(df_exp)}, FDR<0.1: {(df['fdr']<0.1).sum()}")

print("\n=== Enriched (FDR<0.1) ===")
for _, r in df_exp[df_exp["fdr"] < 0.1].iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} log2FC={r['log2FC']:+6.3f} FC={r['median_FC']:.2f} FDR={r['fdr']:.4f} n={r['n_samples']} {s}")

# ── Figure ──
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,
    "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.5,"legend.frameon":False})

df_plot = df_exp.sort_values("log2FC")
fig, ax = plt.subplots(figsize=(4.2, 3.6))
N = len(df_plot)
y = np.arange(N)

ax.axvline(x=0, color="black", linewidth=0.4, zorder=1)

for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr"] < 0.1
    lw = 1.5 if is_sig else 0.7
    alpha = 0.7 if is_sig else 0.25
    ax.plot([0, r["log2FC"]], [i, i], color=color, linewidth=lw, alpha=alpha,
            solid_capstyle="round", zorder=2)

for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "#888888")
    is_sig = r["fdr"] < 0.1
    ax.scatter(r["log2FC"], i, s=42 if is_sig else 24, c=color,
               alpha=0.9 if is_sig else 0.35,
               edgecolors="black" if is_sig else "none", linewidths=0.5 if is_sig else 0, zorder=4)

for i, (_, r) in enumerate(df_plot.iterrows()):
    xi = r["log2FC"]
    if xi >= 0:
        ax.text(xi + 0.06, i, r["gene"], ha="left", va="center", fontsize=6.5, fontstyle="italic", color="#222222")
    else:
        ax.text(xi - 0.06, i, r["gene"], ha="right", va="center", fontsize=6.5, fontstyle="italic", color="#222222")

ax.set_yticks([])
ax.set_xlabel("log2(ILC-high TLS / other TLS)", fontsize=7.5, labelpad=6)
xmax = max(df_plot["log2FC"].max() * 1.3, 0.3)
xmin = min(df_plot["log2FC"].min() * 1.3, -0.3)
ax.set_xlim(xmin, xmax)
ax.tick_params(labelsize=6)

ax.legend([Line2D([0],[0],marker="o",color="w",markerfacecolor="#444444",markersize=7,markeredgecolor="black",markeredgewidth=0.5),
           Line2D([0],[0],marker="o",color="w",markerfacecolor="#bbbbbb",markersize=5,markeredgecolor="none")],
          ["FDR < 0.1", "not significant"], fontsize=6, loc="lower right", handletextpad=0.4)

for j, (cat, color) in enumerate(cat_colors.items()):
    ax.text(0.02, 0.98 - j*0.07, cat, transform=ax.transAxes, fontsize=5.5, color=color, va="top", fontweight="bold")

ax.set_title("ILC-high vs other TLS spots (within-TLS comparison)", fontsize=8, fontweight="bold", loc="left", pad=8)
ax.text(1.0, -0.16, f"n = {n_max} samples | {(df_exp['fdr']<0.1).sum()}/{N} genes FDR < 0.1",
        transform=ax.transAxes, fontsize=5.5, color="#666666", ha="right")

def save_pub(fig, stem, dpi=600):
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.tiff", dpi=dpi, bbox_inches="tight")

fig.tight_layout(pad=0.5)
save_pub(fig, r"E:/GBM/results/fig_receptor_ilc_high_vs_other_tls")
plt.close()
print(f"\nFigure: E:/GBM/results/fig_receptor_ilc_high_vs_other_tls.{{svg,pdf,tiff}}")
