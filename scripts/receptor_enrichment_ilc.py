"""Neurotransmitter receptor enrichment in ILC-dominant TLS spots — category-level analysis"""
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import warnings
warnings.filterwarnings("ignore")

# === Load receptor gene list ===
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
CATEGORIES = {}
cat_names_cn = {}
for col in rx_df.columns:
    name_en = col.strip()
    genes = [g for g in rx_df[col].dropna()]
    CATEGORIES[name_en] = genes
for orig, new in [
    (rx_df.columns[0], "Glutamate"), (rx_df.columns[1], "GABA/Glycine"),
    (rx_df.columns[2], "Cholinergic"), (rx_df.columns[3], "DA/NE"),
    (rx_df.columns[4], "Serotonin"),
]:
    CATEGORIES[new] = CATEGORIES.pop(orig)

ALL_RECEPTORS = sorted(set(g for genes in CATEGORIES.values() for g in genes))
CAT_NAMES = list(CATEGORIES.keys())

TLS_DIR = Path(r"E:\GBM\results\tls_official_cut01")
H5AD_DIR = Path(r"E:\GBM\spatial_data_visium\spatial_data_visium\anndata_with_ilc")
ILC_TYPES = ["ILC1", "ILC2", "ILC3"]
OUT_DIR = Path(r"E:\GBM\results")
OUT_DIR.mkdir(exist_ok=True)

# === Per-sample category scores ===
rows = []
for d in sorted(TLS_DIR.iterdir()):
    if not d.is_dir():
        continue
    tls_csv = d / "tls_spot_scores_official_relaxed.csv"
    h5 = H5AD_DIR / f"{d.name}.h5ad"
    if not tls_csv.exists() or not h5.exists():
        continue

    sid = d.name
    tls = pd.read_csv(tls_csv)
    adata = ad.read_h5ad(h5)
    q05 = adata.obsm["c2l_ilc_q05"]
    if hasattr(q05, "values"):
        q05 = q05.values
    ct = list(adata.uns["c2l_ilc_cell_types"])
    tls_mask = (tls["TLS.region"] == "TLS").values
    dominant = np.array(ct)[q05.argmax(axis=1)]

    # Strict ILC-dominant TLS
    ilc_dom = np.zeros(len(tls), dtype=bool)
    for c in ILC_TYPES:
        idx = ct.index(c)
        p75 = np.percentile(q05[:, idx], 75)
        ilc_dom |= tls_mask & (dominant == c) & (q05[:, idx] >= p75)

    # other TLS = TLS but not ILC-dominant, not ILC-type argmax
    ilc_idx_set = {ct.index(c) for c in ILC_TYPES}
    other_tls = tls_mask & ~ilc_dom & ~np.isin(q05.argmax(axis=1), list(ilc_idx_set))

    if ilc_dom.sum() < 3:
        continue

    # Gene expression
    ge_mask = adata.var["feature_types"] == "Gene Expression"
    ge = adata[:, ge_mask]
    expr = ge.X.toarray() if hasattr(ge.X, "toarray") else ge.X
    var_names = ge.var_names.values

    # per-category score = mean z-score of receptor genes in that category
    for cat, genes in CATEGORIES.items():
        found = [g for g in genes if g in var_names]
        if len(found) < 3:
            continue
        gidx = [list(var_names).index(g) for g in found]
        # raw counts → log1p
        cat_expr = np.log1p(expr[:, gidx])
        # z-score within sample for comparison
        cat_mean = cat_expr.mean(axis=0)
        cat_std = cat_expr.std(axis=0) + 1e-8
        cat_z = ((cat_expr - cat_mean) / cat_std).mean(axis=1)

        for group, mask in [("ILC_TLS", ilc_dom), ("other_TLS", other_tls), ("non_TLS", ~tls_mask)]:
            vals = cat_z[mask]
            if len(vals) == 0:
                continue
            rows.append({"sample": sid, "category": cat, "group": group,
                         "n_spots": len(vals), "score": vals.mean(), "std": vals.std()})

df = pd.DataFrame(rows)
print(f"Samples: {df['sample'].nunique()}")

# === Plot 1: Radar chart — mean category scores per group ===
fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))

group_colors = {"ILC_TLS": "#e41a1c", "other_TLS": "#377eb8", "non_TLS": "#999999"}
angles = np.linspace(0, 2 * np.pi, len(CAT_NAMES), endpoint=False).tolist()
angles += angles[:1]

for group in ["ILC_TLS", "other_TLS", "non_TLS"]:
    means = []
    for cat in CAT_NAMES:
        sub = df[(df["category"] == cat) & (df["group"] == group)]
        if len(sub) > 0:
            means.append(sub["score"].mean())
        else:
            means.append(0)
    means += means[:1]
    ax.fill(angles, means, alpha=0.1, color=group_colors[group])
    ax.plot(angles, means, "o-", linewidth=2, color=group_colors[group], label=group, markersize=6)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(CAT_NAMES, fontsize=11)
ax.set_ylim(-0.1, 0.25)
ax.set_title("Neural receptor category scores by spot type\n(mean log1p z-score)", fontsize=13, pad=25)
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=10)
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_radar_ilc_tls.png", dpi=200, bbox_inches="tight")
plt.close()
print("Radar chart saved")

# === Plot 2: Heatmap — sample × category log2FC (ILC_TLS vs non_TLS) ===
heat_data = {}
sample_order = []
for sid in sorted(df["sample"].unique()):
    sample_order.append(sid)
    for cat in CAT_NAMES:
        ilc = df[(df["sample"] == sid) & (df["category"] == cat) & (df["group"] == "ILC_TLS")]
        non = df[(df["sample"] == sid) & (df["category"] == cat) & (df["group"] == "non_TLS")]
        if len(ilc) > 0 and len(non) > 0:
            fc = ilc["score"].values[0] - non["score"].values[0]
            heat_data[(sid, cat)] = fc

heat_mat = pd.DataFrame(
    [[heat_data.get((s, c), np.nan) for c in CAT_NAMES] for s in sample_order],
    index=sample_order, columns=CAT_NAMES
)
heat_mat = heat_mat.dropna(how="all")

fig, ax = plt.subplots(figsize=(10, max(6, len(heat_mat) * 0.4)))
im = ax.imshow(heat_mat.values, cmap="RdBu_r", aspect="auto", vmin=-0.3, vmax=0.3)
ax.set_xticks(range(len(CAT_NAMES)))
ax.set_xticklabels(CAT_NAMES, fontsize=10, rotation=45, ha="right")
ax.set_yticks(range(len(heat_mat)))
ax.set_yticklabels(heat_mat.index, fontsize=8)
ax.set_title("ILC-TLS vs non-TLS: receptor category score difference\n(Δ z-score, red = enriched)", fontsize=12)
plt.colorbar(im, ax=ax, shrink=0.8)
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_heatmap_ilc_tls.png", dpi=200, bbox_inches="tight")
plt.close()
print("Heatmap saved")

# === Plot 3: Boxplot — category-level comparison ===
fig, axes = plt.subplots(1, 5, figsize=(18, 5), sharey=True)
for i, (cat, ax) in enumerate(zip(CAT_NAMES, axes)):
    data_groups = []
    labels = []
    for group, color in [("ILC_TLS", "#e41a1c"), ("other_TLS", "#377eb8"), ("non_TLS", "#999999")]:
        vals = df[(df["category"] == cat) & (df["group"] == group)]["score"].dropna().values
        if len(vals) > 0:
            data_groups.append(vals)
            labels.append(group)
    bp = ax.boxplot(data_groups, labels=labels, patch_artist=True, widths=0.6)
    for patch, color in zip(bp["boxes"], [group_colors[l] for l in labels]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
    ax.set_title(cat, fontsize=11)
    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)

    # Add p-values
    if len(data_groups) >= 2:
        _, p = mannwhitneyu(data_groups[0], data_groups[2], alternative="two-sided")
        n_ilc = len(data_groups[0])
        ax.text(1, ax.get_ylim()[1] * 0.95, f"p={p:.3f}\nn={n_ilc}", ha="center", fontsize=8)

axes[0].set_ylabel("Receptor category z-score", fontsize=11)
fig.suptitle("Neurotransmitter receptor expression by spot type", fontsize=13)
fig.tight_layout()
fig.savefig(OUT_DIR / "receptor_boxplot_ilc_tls.png", dpi=200, bbox_inches="tight")
plt.close()
print("Boxplot saved")

# === Summary stats ===
print("\n=== Category-level enrichment ===")
for cat in CAT_NAMES:
    ilc = df[(df["category"] == cat) & (df["group"] == "ILC_TLS")]["score"]
    non = df[(df["category"] == cat) & (df["group"] == "non_TLS")]["score"]
    if len(ilc) >= 3 and len(non) >= 3:
        diff = ilc.mean() - non.mean()
        _, p = mannwhitneyu(ilc, non, alternative="two-sided")
        print(f"  {cat:15s}: Δz={diff:+.4f}  p={p:.4f}  n_samples={len(ilc)}")

from scipy.stats import false_discovery_control
pvals = []
for cat in CAT_NAMES:
    ilc = df[(df["category"] == cat) & (df["group"] == "ILC_TLS")]["score"]
    non = df[(df["category"] == cat) & (df["group"] == "non_TLS")]["score"]
    if len(ilc) >= 3 and len(non) >= 3:
        _, p = mannwhitneyu(ilc, non, alternative="two-sided")
        pvals.append((cat, p))
if pvals:
    fdr_vals = false_discovery_control([p for _, p in pvals])
    print("\nFDR-corrected:")
    for (cat, _), fdr in zip(pvals, fdr_vals):
        print(f"  {cat:15s}: FDR={fdr:.4f}")

print(f"\nSaved to {OUT_DIR}")
