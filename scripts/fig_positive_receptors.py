"""
Positive receptor definition + 4-panel Nature figure
Criteria: n>=20, pct_expr>=5%, delta>0, FDR<0.1, prop(delta>0)>=60%
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
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

# ====== Compute per-sample per-gene: delta, log2FC, pct_expr ======
gene_sample_data = {}  # gene -> [{sample, delta, log2fc, pct_ie, pct_ne, pct_ot}]
total_samples = 0

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

        # ILC-enriched TLS
        ilc_high = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx[c]] >= GLOBAL_THRESH[c])
        ilc_enr = ilc_high & tls_mask
        other_tls = tls_mask & ~ilc_enr
        if ilc_enr.sum() < 3: continue

        # Gene expression
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_idx = [list(vn_all).index(g) for g in rx_present]
        rx_dense = ge[:, rx_idx].X.toarray() if hasattr(ge[:, rx_idx].X,"toarray") else np.asarray(ge[:, rx_idx].X)
        expr_norm = np.log1p(rx_dense / (lib_size[:,None]/10000))

        ie_expr = expr_norm[ilc_enr]; ne_expr = expr_norm[~tls_mask]; ot_expr = expr_norm[other_tls]
        for gi, g in enumerate(rx_present):
            ie_m = ie_expr[:, gi].mean(); ne_m = ne_expr[:, gi].mean(); ot_m = ot_expr[:, gi].mean()
            if ie_m == 0 or ne_m == 0: continue
            delta = ie_m - ne_m
            l2fc = np.log2(ie_m / ne_m)
            pct_ie = (rx_dense[ilc_enr, gi] > 0).mean()
            pct_ne = (rx_dense[~tls_mask, gi] > 0).mean()
            pct_ot = (rx_dense[other_tls, gi] > 0).mean()
            gene_sample_data.setdefault(g, []).append({"sample":sd.name, "delta":delta, "log2fc":l2fc,
                "pct_ie":pct_ie, "pct_ne":pct_ne, "pct_ot":pct_ot})
        total_samples += 1

print(f"Total samples: {total_samples}")

# ====== Aggregate per gene with multi-criteria ======
rows = []
for g, vals in gene_sample_data.items():
    n = len(vals); deltas = np.array([v["delta"] for v in vals])
    l2fcs = np.array([v["log2fc"] for v in vals])
    pct_ies = np.array([v["pct_ie"] for v in vals])
    pct_ots = np.array([v["pct_ot"] for v in vals])
    pct_nes = np.array([v["pct_ne"] for v in vals])

    median_delta = np.median(deltas); median_l2fc = np.median(l2fcs)
    median_pct_ie = np.median(pct_ies); median_pct_ot = np.median(pct_ots)
    prop_pos = (deltas > 0).mean()

    _, p_delta = wilcoxon(deltas, alternative="two-sided")
    _, p_l2fc = wilcoxon(l2fcs, alternative="two-sided")

    rows.append({"gene":g, "n":n, "median_delta":median_delta, "median_log2FC":median_l2fc,
        "pct_ILC_TLS":median_pct_ie, "pct_other_TLS":median_pct_ot, "pct_nonTLS":pct_nes[0],
        "prop_positive":prop_pos, "p_delta":p_delta, "p_log2fc":p_l2fc,
        "category":gene_cat.get(g,"?")})

df_all = pd.DataFrame(rows)
df_all["fdr_delta"] = false_discovery_control(df_all["p_delta"].values)
df_all["fdr_log2fc"] = false_discovery_control(df_all["p_log2fc"].values)
# Fill NaN pct values before use in plots
df_all["pct_other_TLS"] = df_all["pct_other_TLS"].fillna(0)
df_all["pct_nonTLS"] = df_all["pct_nonTLS"].fillna(0)

# ====== Positive receptor criteria ======
POSITIVE = (
    (df_all["n"] >= 20) &
    (df_all["pct_ILC_TLS"] >= 0.05) &
    (df_all["median_delta"] > 0) &
    (df_all["fdr_delta"] < 0.1) &
    (df_all["prop_positive"] >= 0.6)
)
df_all["positive"] = POSITIVE
df_pos = df_all[POSITIVE].sort_values("median_delta", ascending=False)
df_all.to_csv(OUT / "receptor_multi_criteria.csv", index=False)

print(f"\nPositive receptors (5 criteria): {len(df_pos)}")
for _, r in df_pos.iterrows():
    print(f"  {r['gene']:8s} {r['category']:22s} delta={r['median_delta']:+.4f} pct={r['pct_ILC_TLS']:.3f} prop+={r['prop_positive']:.2f} n={r['n']}")

# ====== Figure: 4-panel ======
fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(10, 8))

# ---- Panel A: Volcano / ranking ----
# (ax_a already created)
df_a = df_all.copy()
for _, r in df_a.iterrows():
    c = CAT.get(r["category"], "#888")
    sz = 15 + max(r["pct_ILC_TLS"], 0.01) * 120
    is_pos = r["positive"]
    ax_a.scatter(r["median_delta"], -np.log10(max(r["p_delta"], 1e-10)), s=sz, c=c,
                 alpha=0.85 if is_pos else 0.3, edgecolors="black" if is_pos else "none", linewidths=0.4, zorder=5 if is_pos else 2)
    if is_pos:
        ax_a.annotate(r["gene"], (r["median_delta"], -np.log10(max(r["p_delta"], 1e-10))),
                      fontsize=5.5, fontstyle="italic", xytext=(3, 3), textcoords="offset points")
ax_a.axhline(y=-np.log10(0.1), color="#ccc", linewidth=0.5, linestyle="--")
ax_a.axvline(x=0, color="#ccc", linewidth=0.5)
ax_a.set_xlabel("median delta (ILC-enriched TLS - non-TLS)", fontsize=7)
ax_a.set_ylabel("-log10(p)", fontsize=7)
ax_a.set_title(f"a  Volcano: {len(df_pos)} positive receptors\n(5 criteria: n>=20, pct>=5%, delta>0, FDR<0.1, prop+>=60%)",
               fontsize=8, fontweight="bold", loc="left")
size_legend = [Line2D([0],[0], marker='o', color='w', markerfacecolor='gray', markersize=np.sqrt(s)/2, alpha=0.4)
               for s in [25, 55, 85]]
ax_a.legend(size_legend, ["10%", "30%", "50%"], fontsize=5.5, loc="upper right", title="ILC-TLS rate", title_fontsize=6, borderpad=0.3)

# ---- Panel B: Dot plot: receptor x group ----
df_b = df_pos.nlargest(15, "median_delta")
groups = ["nonTLS", "other TLS", "ILC-TLS"]
y_ticks = np.arange(len(df_b))
for gi, (_, r) in enumerate(df_b.iterrows()):
    vals = [r["pct_nonTLS"], r.get("pct_other_TLS", 0) or 0, r["pct_ILC_TLS"]]
    sizes = [15 + max(v, 0.01) * 150 for v in vals]
    colors_dot = ["#ccc", "#4c72b0", "#c44e52"]
    for gj in range(3):
        ax_b.scatter(gj, gi, s=sizes[gj], c=colors_dot[gj], alpha=0.8, edgecolors="black", linewidths=0.3, zorder=5)
ax_b.set_yticks(y_ticks)
ax_b.set_yticklabels(df_b["gene"].values, fontsize=6.5, fontstyle="italic")
ax_b.set_xticks(range(3))
ax_b.set_xticklabels(groups, fontsize=6)
ax_b.set_xlim(-0.5, 2.5)
ax_b.set_title(f"b  Detection rate by group (top {len(df_b)} positive)", fontsize=8, fontweight="bold", loc="left")
ax_b.tick_params(axis='x', rotation=30)

# ---- Panel C: Paired sample-level for top 4 ----
top4 = [g for g in df_pos.head(4)["gene"].values if g in gene_sample_data and len(gene_sample_data[g]) > 5]
for gi, g in enumerate(top4[:4]):  # max 4
    vals = gene_sample_data[g]
    ie_vals = np.array([v["pct_ie"] for v in vals])
    ot_vals = np.array([v["pct_ot"] for v in vals])
    # Plot paired lines (thin, transparent)
    for i in range(min(len(ie_vals), 100)):  # limit to 100 lines for performance
        ax_c.plot([gi - 0.25, gi + 0.25], [ot_vals[i], ie_vals[i]],
                  color="#ccc", linewidth=0.2, alpha=0.3, zorder=2)
    # Individual points
    ax_c.scatter(np.full(len(ot_vals), gi - 0.25), ot_vals, s=6, c="#4c72b0", alpha=0.3, zorder=3)
    ax_c.scatter(np.full(len(ie_vals), gi + 0.25), ie_vals, s=6, c="#c44e52", alpha=0.3, zorder=3)
    # Medians
    ax_c.scatter(gi - 0.25, np.median(ot_vals), s=40, c="#4c72b0", edgecolors="black", linewidths=0.5, zorder=6)
    ax_c.scatter(gi + 0.25, np.median(ie_vals), s=40, c="#c44e52", edgecolors="black", linewidths=0.5, zorder=6)
ax_c.set_xticks(range(len(top4[:4])))
ax_c.set_xticklabels(top4[:4], fontsize=7, fontstyle="italic")
ax_c.set_ylabel("pct_expr (ILC-enriched TLS)", fontsize=7)
ax_c.set_title("c  Paired per-sample: other TLS vs ILC-enriched TLS", fontsize=8, fontweight="bold", loc="left")
ax_c.legend([Line2D([0],[0],marker='o',c='w',markerfacecolor='#4c72b0', markersize=8),
             Line2D([0],[0],marker='o',c='w',markerfacecolor='#c44e52', markersize=8)],
            ["other TLS", "ILC-TLS"], fontsize=6, loc="upper left")

# ---- Panel D: Category summary ----
cat_counts = df_pos["category"].value_counts()
cat_colors_plot = [CAT.get(c, "#888") for c in cat_counts.index]
bars = ax_d.barh(range(len(cat_counts)), cat_counts.values, color=cat_colors_plot, alpha=0.6, edgecolor="black", linewidth=0.3)
# Also show non-positive count per category
for i, (cat, n_pos) in enumerate(cat_counts.items()):
    n_total = sum(1 for _, r in df_all.iterrows() if r["category"] == cat)
    ax_d.text(n_pos + 0.3, i, f"{n_pos}/{n_total}", fontsize=6, va="center", color="#666")
ax_d.set_yticks(range(len(cat_counts)))
ax_d.set_yticklabels([c.split("(")[0].strip() for c in cat_counts.index], fontsize=6)
ax_d.set_xlabel("positive receptors", fontsize=7)
ax_d.set_title(f"d  Positive receptors by category ({len(df_pos)} total)", fontsize=8, fontweight="bold", loc="left")

fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_positive_receptors_4panel")
plt.close()
print(f"\nFigure: {OUT}/fig_positive_receptors_4panel.{{svg,pdf,tiff}}")
