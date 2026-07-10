"""Per-gene neurotransmitter receptor enrichment in ILC-dominant TLS spots"""
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
CATEGORIES = {}
for col in rx_df.columns:
    genes = [g for g in rx_df[col].dropna()]
    for g in genes:
        CATEGORIES[g] = col.strip()

ALL_RECEPTORS = sorted(CATEGORIES.keys())

TLS_DIR = Path(r"E:\GBM\results\tls_official_cut01")
H5AD_DIR = Path(r"E:\GBM\spatial_data_visium\spatial_data_visium\anndata_with_ilc")
ILC_TYPES = ["ILC1", "ILC2", "ILC3"]
OUT_DIR = Path(r"E:\GBM\results")

# Per-sample, per-gene: mean log1p expression in ILC-TLS vs non-TLS
gene_data = {}
ilc_counts = {}

for d in sorted(TLS_DIR.iterdir()):
    if not d.is_dir():
        continue
    tls_csv = d / "tls_spot_scores_official_relaxed.csv"
    h5 = H5AD_DIR / f"{d.name}.h5ad"
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

    ge_mask = adata.var["feature_types"] == "Gene Expression"
    ge = adata[:, ge_mask]
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
            ilc_counts[gene] = []
        gene_data[gene].append(fc)
        ilc_counts[gene].append(ilc_dom.sum())

# Aggregate per gene
rows = []
for gene, fcs in gene_data.items():
    if len(fcs) < 5:
        continue
    arr = np.array(fcs)
    arr = arr[~np.isnan(arr) & ~np.isinf(arr)]
    if len(arr) < 5:
        continue
    median_fc = np.median(arr)
    mean_fc = np.mean(arr)
    # paired test: FC differs from 1.0
    _, pval = wilcoxon(arr - 1.0, alternative="two-sided")
    rows.append({
        "gene": gene, "n_samples": len(arr),
        "median_FC": median_fc, "mean_FC": mean_fc,
        "log2FC": np.log2(median_fc), "pvalue": pval,
        "total_ilc_spots": sum(ilc_counts[gene]),
        "category": CATEGORIES.get(gene, "?"),
    })

df = pd.DataFrame(rows)
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("log2FC", ascending=False)

# Print significant genes
sig = df[df["fdr"] < 0.1]
print(f"Total genes tested: {len(df)}")
print(f"FDR < 0.1: {len(sig)}")
print(f"\nSignificant genes (FDR < 0.1):")
for _, r in sig.iterrows():
    stars = "***" if r["fdr"] < 0.001 else "**" if r["fdr"] < 0.01 else "*"
    print(f"  {r['gene']:10s}  log2FC={r['log2FC']:+.3f}  FC={r['median_FC']:.2f}  FDR={r['fdr']:.4f}  n={r['n_samples']}  {stars}")

# Volcano plot
fig, ax = plt.subplots(figsize=(12, 8))
colors = {"兴奋性受体": "#E64A19", "抑制性受体": "#1B5E20", "胆碱能受体": "#0D47A1",
          "多巴胺/去甲肾上腺素": "#6A1B9A", "5-羟色胺受体": "#BF360C"}
cat_map = {rx_df.columns[0]: "兴奋性受体", rx_df.columns[1]: "抑制性受体",
           rx_df.columns[2]: "胆碱能受体", rx_df.columns[3]: "多巴胺/去甲肾上腺素",
           rx_df.columns[4]: "5-羟色胺受体"}

for _, r in df.iterrows():
    cat_cn = cat_map.get(r["category"], r["category"])
    color = colors.get(cat_cn, "gray")
    alpha = 0.8 if r["fdr"] < 0.1 else 0.3
    size = 40 if r["fdr"] < 0.05 else 20
    ax.scatter(r["log2FC"], -np.log10(r["pvalue"]), c=color, alpha=alpha, s=size,
               edgecolors="black" if r["fdr"] < 0.1 else "none", linewidths=0.5)

# Label significant + notable genes
for _, r in df.iterrows():
    if r["fdr"] < 0.1 or abs(r["log2FC"]) > 0.5:
        ax.annotate(r["gene"], (r["log2FC"], -np.log10(r["pvalue"])),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points")

ax.axhline(y=-np.log10(0.1), color="gray", linestyle="--", alpha=0.5, label="FDR=0.1")
ax.axvline(x=0, color="gray", linestyle="--", alpha=0.3)
ax.set_xlabel("log2 Fold Change (ILC-TLS / non-TLS)", fontsize=12)
ax.set_ylabel("-log10(p-value)", fontsize=12)
ax.set_title("Neurotransmitter receptors: ILC-dominant TLS vs non-TLS\n", fontsize=14)

# Legend
from matplotlib.patches import Patch
legend_elems = [Patch(facecolor=colors[k], alpha=0.6, label=k) for k in colors]
ax.legend(handles=legend_elems, fontsize=9, loc="upper right")
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_volcano_ilc_tls.png", dpi=200, bbox_inches="tight")
plt.close()

# Dot plot (top genes)
top = df.nlargest(20, "log2FC")
fig, ax = plt.subplots(figsize=(10, 6))
ypos = range(len(top))
ax.scatter(top["log2FC"].values, ypos, c=[colors.get(cat_map.get(c, c), "gray") for c in top["category"]],
           s=top["n_samples"] * 8, alpha=0.8, edgecolors="black", linewidths=0.5)
ax.set_yticks(ypos)
ax.set_yticklabels(top["gene"].values, fontsize=10)
ax.axvline(x=0, color="gray", linestyle="--")
ax.set_xlabel("log2FC (ILC-TLS / non-TLS)", fontsize=12)
ax.set_title("Top 20 enriched neurotransmitter receptors\nin ILC-dominant TLS spots", fontsize=13)
# Size legend
for ns in [5, 10, 15]:
    ax.scatter([], [], s=ns * 8, c="gray", alpha=0.5, label=f"n={ns}")
ax.legend(title="samples", fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_dotplot_ilc_tls.png", dpi=200, bbox_inches="tight")
plt.close()

print(f"\nSaved to {OUT_DIR}")
df.to_csv(OUT_DIR / "receptor_per_gene_ilc.csv", index=False)
