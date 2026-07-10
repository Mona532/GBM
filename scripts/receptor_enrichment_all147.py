"""Per-gene receptor enrichment across all 147 samples (dataset 1 + dataset 2)"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

# === Config ===
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)", "Inhib (GABA/Gly)", "Cholinergic (ACh)", "DA/NE", "Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna():
        gene_cat[g] = cat_order[idx]
cat_colors = {"Excit (Glutamate)": "#C24B3C", "Inhib (GABA/Gly)": "#358554",
              "Cholinergic (ACh)": "#3575A3", "DA/NE": "#764E9F", "Serotonin (5-HT)": "#BF5A2E"}

ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1", "ILC2", "ILC3"]

# === Dataset definitions ===
DATASETS = [
    {
        "name": "dryad_h5ad",
        "tls_dir": Path(r"E:/GBM/results/tls_official_cut01"),
        "h5ad_dir": Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc"),
    },
    {
        "name": "visium_all",
        "tls_dir": Path(r"E:/GBM/results/tls_visium_all"),
        "h5ad_dir": Path(r"E:/GBM/ST_DATA/visium_all_h5ad"),
    },
]
# Exclude DMG samples (pediatric, not GBM)
SKIP_DMG = True

OUT = Path(r"E:/GBM/results")
OUT.mkdir(exist_ok=True)

# === Per-sample analysis ===
gene_data = {}

for ds in DATASETS:
    for d in sorted(ds["tls_dir"].iterdir()):
        if not d.is_dir():
            continue
        tls_csv = d / "tls_spot_scores_official_relaxed.csv"
        h5 = ds["h5ad_dir"] / f"{d.name}.h5ad"
        if not tls_csv.exists() or not h5.exists():
            continue
        tls = pd.read_csv(tls_csv)
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05, "values"):
            q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        tls_mask = (tls["TLS.region"] == "TLS").values
        dominant = np.array(ct)[q05.argmax(axis=1)]

        ilc_dom = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES:
            idx = ct.index(c)
            p75 = np.percentile(q05[:, idx], 75)
            ilc_dom |= tls_mask & (dominant == c) & (q05[:, idx] >= p75)

        if ilc_dom.sum() < 3:
            continue

        # Gene expression - handle both formats
        if "feature_types" in adata.var:
            ge_mask = adata.var["feature_types"] == "Gene Expression"
            ge = adata[:, ge_mask]
        else:
            ge = adata
        expr = ge.X.toarray() if hasattr(ge.X, "toarray") else ge.X
        var_names = ge.var_names.values

        for gene in ALL_RECEPTORS:
            if gene not in var_names:
                continue
            gidx = list(var_names).index(gene)
            ilc_mean = np.log1p(expr[ilc_dom, gidx]).mean()
            ctrl_mean = np.log1p(expr[~tls_mask, gidx]).mean()
            fc = ilc_mean / ctrl_mean if ctrl_mean > 0 else np.nan
            if gene not in gene_data:
                gene_data[gene] = []
            gene_data[gene].append(fc)

# === Aggregate per gene ===
rows = []
for gene, fcs in gene_data.items():
    arr = np.array([x for x in fcs if not np.isnan(x) and not np.isinf(x)])
    if len(arr) < 5:
        continue
    _, pval = wilcoxon(arr - 1.0, alternative="two-sided")
    rows.append({"gene": gene, "n_samples": len(arr),
                 "median_FC": np.median(arr), "log2FC": np.log2(np.median(arr)),
                 "pvalue": pval, "category": gene_cat.get(gene, "?")})

df = pd.DataFrame(rows)
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("log2FC", ascending=False)
df.to_csv(OUT / "receptor_all147_per_gene.csv", index=False)

# Split: expressed (log2FC > -10) vs not expressed in ILC-TLS (FC=0)
df_exp = df[df["log2FC"] > -10].copy()
df_zero = df[df["log2FC"] <= -10].copy()
n_max = df["n_samples"].max()
print(f"Total genes: {len(df)} | Expressed in ILC-TLS: {len(df_exp)} | NOT expressed in ILC-TLS: {len(df_zero)}")
print(f"Samples with ILC-dominant TLS: {n_max}\n")

print("=== ENRICHED in ILC-TLS (FDR < 0.1) ===")
for _, r in df_exp[df_exp["fdr"] < 0.1].iterrows():
    s = "***" if r["fdr"] < 0.001 else "**" if r["fdr"] < 0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} log2FC={r['log2FC']:+7.3f}  FC={r['median_FC']:.2f}  FDR={r['fdr']:.4f}  n={r['n_samples']}  {s}")

print("\n=== DEPLETED / NOT EXPRESSED in ILC-TLS (ILC_mean=0, FDR < 0.1) ===")
df_zero_sig = df_zero[df_zero["fdr"] < 0.1].sort_values("n_samples", ascending=False)
for _, r in df_zero_sig.iterrows():
    s = "***" if r["fdr"] < 0.001 else "**" if r["fdr"] < 0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} FC=0  FDR={r['fdr']:.4f}  n={r['n_samples']}  {s}")

# === Figure: two-panel (A: expressed lollipop, B: zero-expression summary) ===
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.5,
    "legend.frameon": False,
})

fig = plt.figure(figsize=(8, 5))

# Panel A: Expressed genes lollipop
ax_a = fig.add_axes([0.08, 0.12, 0.52, 0.82])
df_plot = df_exp.sort_values("log2FC", ascending=True)
y = np.arange(len(df_plot))
x = df_plot["log2FC"].values
colors_plot = [cat_colors.get(r["category"], "gray") for _, r in df_plot.iterrows()]

ax_a.axvline(x=0, color="black", linewidth=0.6, zorder=2)
for i, (_, r) in enumerate(df_plot.iterrows()):
    ci = cat_colors.get(r["category"], "gray")
    alpha = 0.15 if r["fdr"] >= 0.1 else 0.6
    ax_a.plot([0, r["log2FC"]], [i, i], color=ci, linewidth=1.2 if alpha > 0.3 else 0.5,
              alpha=alpha, zorder=3, solid_capstyle="round")
for i, (_, r) in enumerate(df_plot.iterrows()):
    color = cat_colors.get(r["category"], "gray")
    is_sig = r["fdr"] < 0.1
    ax_a.scatter(r["log2FC"], i, s=32 if is_sig else 16, c=color, alpha=0.9 if is_sig else 0.3,
                 edgecolors="black" if is_sig else "none", linewidths=0.4, zorder=5)
margin = abs(x.min()) * 0.15
for i, gene in enumerate(df_plot["gene"]):
    ax_a.text(x.min() - margin, i, gene, ha="right", va="center", fontsize=5.5, fontstyle="italic")
ax_a.set_yticks([])
ax_a.set_xlabel("log2(ILC-TLS / non-TLS)", fontsize=7)
ax_a.set_title(f"a  Expressed receptors (n={n_max} samples)", fontsize=8, fontweight="bold", loc="left")
ax_a.tick_params(labelsize=5.5)
ax_a.spines["left"].set_visible(False)
from matplotlib.lines import Line2D
ax_a.legend([Line2D([0],[0],marker="o",color="w",markerfacecolor="gray",markersize=6,markeredgecolor="black",markeredgewidth=0.4),
             Line2D([0],[0],marker="o",color="w",markerfacecolor="gray",markersize=3.5)],
            ["FDR<0.1","NS"], fontsize=5, loc="lower right", handletextpad=0.5)

# Panel B: Zero-expression genes by category
ax_b = fig.add_axes([0.64, 0.12, 0.34, 0.82])
df_z = df_zero.copy()
# Count zero genes per category
cat_counts = {}
for cat in cat_order:
    n_total = sum(1 for g in gene_cat if gene_cat[g] == cat)
    n_exp = sum(1 for _, r in df_exp.iterrows() if r["category"] == cat)
    n_zero = sum(1 for _, r in df_zero.iterrows() if r["category"] == cat)
    n_sig_zero = sum(1 for _, r in df_zero_sig.iterrows() if r["category"] == cat)
    cat_counts[cat] = (n_total, n_exp, n_sig_zero, n_zero)

cats = list(cat_counts.keys())
n_tot = [cat_counts[c][0] for c in cats]
n_exp = [cat_counts[c][1] for c in cats]
n_sig_zero = [cat_counts[c][2] for c in cats]
n_ns_zero = [cat_counts[c][3] - cat_counts[c][2] for c in cats]
colors_cat = [cat_colors[c] for c in cats]

y_pos = range(len(cats))
bar_exp = ax_b.barh(y_pos, n_exp, height=0.5, color=colors_cat, alpha=0.7, label="Expressed in ILC-TLS")
bar_sig = ax_b.barh(y_pos, n_sig_zero, height=0.5, left=n_exp, color=colors_cat, alpha=0.4,
                     edgecolor="black", linewidth=0.5, hatch="///", label="Zero (FDR<0.1)")
bar_ns = ax_b.barh(y_pos, n_ns_zero, height=0.5, left=np.array(n_exp)+np.array(n_sig_zero),
                    color=colors_cat, alpha=0.12, label="Zero (NS)")

# Labels
for i, cat in enumerate(cats):
    total = n_tot[i]
    ax_b.text(total + 1, i, f"{total}", fontsize=6, va="center", fontweight="bold")

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(cats, fontsize=6)
ax_b.set_xlabel("Number of genes", fontsize=7)
ax_b.set_title(f"b  Not expressed in ILC-TLS", fontsize=8, fontweight="bold", loc="left")
ax_b.legend(fontsize=5.5, loc="lower right", handletextpad=0.5)
ax_b.set_xlim(0, max(n_tot) * 1.3)

fig.savefig(OUT / "fig_receptor_all147.svg", bbox_inches="tight", dpi=600)
fig.savefig(OUT / "fig_receptor_all147.pdf", bbox_inches="tight", dpi=600)
fig.savefig(OUT / "fig_receptor_all147.tiff", bbox_inches="tight", dpi=600)
plt.close()
print(f"\nSaved: fig_receptor_all147.{{svg,pdf,tiff}}")
print(f"CSV: {OUT}/receptor_all147_per_gene.csv")
