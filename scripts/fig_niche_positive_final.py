"""
Niche-stratified positive receptor identification (CellCharter global niches)
ILC-enriched TLS vs other TLS within same niche, then aggregate across niches
5+1 criteria: n>=20, pct>=5%, delta>0, FDR<0.1, prop+>=60%, pct_ILC > pct_other
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import gc, warnings; warnings.filterwarnings("ignore")

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

# ====== Step 1: Global CellCharter niche ======
adata_list = []
sample_records = []  # (sid, n_spots, tls_mask, expr_norm, vn_rx, q05, ct_names, ilc_idx, coords, rx_idx_map)
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

        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_idx_list = [list(vn_all).index(g) for g in rx_present]
        rx_raw = ge[:, rx_idx_list].X.toarray() if hasattr(ge[:, rx_idx_list].X,"toarray") else np.asarray(ge[:, rx_idx_list].X)
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))

        # Fix 1: per-sample z-scored q05 for CellCharter
        q05_z = (q05 - q05.mean(axis=0)) / (q05.std(axis=0) + 1e-8)
        obs_df = pd.DataFrame(index=adata.obs_names); obs_df["sample"] = sd.name
        a = ad.AnnData(X=q05_z.astype(np.float32), obs=obs_df, obsm={"spatial": coords})
        a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
        adata_list.append(a)
        sample_records.append((sd.name, q05.shape[0], tls_mask, expr_norm, rx_raw, rx_present, q05, ct_names, ilc_idx_map, coords))

adata_all = ad.concat(adata_list, join="inner", label="sample_id")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adata_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
autok = cc.tl.ClusterAutoK(n_clusters=(5, 20), max_runs=5)
autok.fit(adata_all, use_rep="X_cellcharter")
niche_labels_all = autok.predict(adata_all, use_rep="X_cellcharter")
K = len(set(niche_labels_all))
print(f"Global niches: {K}")

# ====== Step 2: Per-sample × niche per-gene comparison ======
gene_niche_data = {}  # gene -> [{sample, niche, delta, pct_ie, pct_ot, l2fc}]
offset = 0
for sid, n_spots, tls_mask, expr_norm, rx_raw, rx_present, q05, ct_names, ilc_idx_map, coords in sample_records:
    niche_labels = niche_labels_all[offset:offset+n_spots]
    offset += n_spots

    for niche_id in sorted(set(niche_labels)):
        niche_mask = niche_labels == niche_id
        tls_in_niche = tls_mask & niche_mask
        if tls_in_niche.sum() < 10: continue

        ilc_enr = np.zeros(tls_in_niche.sum(), dtype=bool)
        tls_spot_idx = np.where(tls_in_niche)[0]
        for c in ILC_TYPES:
            ilc_enr |= (q05[tls_spot_idx, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
        other_tls = ~ilc_enr

        if ilc_enr.sum() < 3 or other_tls.sum() < 3: continue

        ie_expr = expr_norm[tls_spot_idx[ilc_enr]]
        ot_expr = expr_norm[tls_spot_idx[other_tls]]

        for gi, g in enumerate(rx_present):
            ie_m = ie_expr[:, gi].mean(); ot_m = ot_expr[:, gi].mean()
            delta = ie_m - ot_m
            l2fc = np.log2((ie_m + 1e-6) / (ot_m + 1e-6))
            # pct_expr from raw counts > 0
            pct_ie = (rx_raw[tls_spot_idx[ilc_enr], gi] > 0).mean()
            pct_ot = (rx_raw[tls_spot_idx[other_tls], gi] > 0).mean()
            gene_niche_data.setdefault(g, []).append({
                "sample": sid, "niche": niche_id, "delta": delta, "log2fc": l2fc,
                "pct_ie": pct_ie, "pct_ot": pct_ot})

print(f"Total niche-gene pairs: {sum(len(v) for v in gene_niche_data.values())}")

# ====== Step 3: Aggregate per gene, per-sample (median across niches) ======
gene_sample = {}
for g, niche_vals in gene_niche_data.items():
    # Group by sample, compute median across niches
    sample_vals = {}
    for v in niche_vals:
        s = v["sample"]
        sample_vals.setdefault(s, []).append(v)
    for s, vals in sample_vals.items():
        gene_sample.setdefault(g, []).append({
            "sample": s, "delta": np.median([v["delta"] for v in vals]),
            "log2fc": np.median([v["log2fc"] for v in vals]),
            "pct_ie": np.median([v["pct_ie"] for v in vals]),
            "pct_ot": np.median([v["pct_ot"] for v in vals])})

rows = []
for g, svals in gene_sample.items():
    n = len(svals); deltas = np.array([v["delta"] for v in svals])
    pct_ies = np.array([v["pct_ie"] for v in svals])
    pct_ots = np.array([v["pct_ot"] for v in svals])
    l2fcs = np.array([v["log2fc"] for v in svals])

    med_delta = np.median(deltas); med_l2fc = np.median(l2fcs)
    med_pct_ie = np.median(pct_ies); med_pct_ot = np.median(pct_ots)
    prop_pos = (deltas > 0).mean()

    _, p_delta = wilcoxon(deltas, alternative="two-sided")
    _, p_l2fc = wilcoxon(l2fcs, alternative="two-sided")

    rows.append({"gene":g, "n":n, "median_delta":med_delta, "median_log2FC":med_l2fc,
        "pct_ILC_TLS":med_pct_ie, "pct_other_TLS":med_pct_ot,
        "prop_positive":prop_pos, "p_delta":p_delta, "p_l2fc":p_l2fc,
        "category":gene_cat.get(g,"?")})

df_all = pd.DataFrame(rows)
df_all = df_all[df_all["p_delta"].notna() & (df_all["p_delta"] > 0) & (df_all["p_delta"] < 1)]
df_all["fdr_delta"] = false_discovery_control(df_all["p_delta"].values)
df_all["fdr_log2fc"] = false_discovery_control(df_all["p_l2fc"].values)

# ====== Step 4: Positive receptor criteria ======
POSITIVE = (
    (df_all["n"] >= 20) &
    (df_all["pct_ILC_TLS"] >= 0.05) &
    (df_all["median_delta"] > 0) &
    (df_all["fdr_delta"] < 0.1) &
    (df_all["prop_positive"] >= 0.6) &
    (df_all["pct_ILC_TLS"] > df_all["pct_other_TLS"])
)
df_all["positive"] = POSITIVE
df_pos = df_all[POSITIVE].sort_values("median_delta", ascending=False)
df_all.to_csv(OUT / "receptor_niche_positive.csv", index=False)

print(f"Positive receptors (niche-stratified, 6 criteria): {len(df_pos)}")
for _, r in df_pos.iterrows():
    print(f"  {r['gene']:8s} {r['category']:22s} delta={r['median_delta']:+.4f} pct_IE={r['pct_ILC_TLS']:.3f} pct_OT={r['pct_other_TLS']:.3f} prop+={r['prop_positive']:.2f} n={r['n']}")

# ====== Figure: 4-panel ======
fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(10, 8))

# Panel A: Volcano
for _, r in df_all.iterrows():
    c = CAT.get(r["category"], "#888")
    sz = 15 + max(r["pct_ILC_TLS"], 0.01) * 120
    is_pos = r["positive"]
    ax_a.scatter(r["median_delta"], -np.log10(max(r["p_delta"], 1e-10)), s=sz, c=c,
                 alpha=0.85 if is_pos else 0.3, edgecolors="black" if is_pos else "none", linewidths=0.4, zorder=5 if is_pos else 2)
    if is_pos:
        ax_a.annotate(r["gene"], (r["median_delta"], -np.log10(max(r["p_delta"], 1e-10))),
                      fontsize=5, fontstyle="italic", xytext=(3, 3), textcoords="offset points")
ax_a.axhline(y=-np.log10(0.1), color="#ccc", linewidth=0.5, linestyle="--")
ax_a.axvline(x=0, color="#ccc", linewidth=0.5)
ax_a.set_xlabel("median delta (ILC-enriched - other TLS)", fontsize=7)
ax_a.set_ylabel("-log10(p)", fontsize=7)
ax_a.set_title(f"a  Niche-stratified volcano ({len(df_pos)} positive)\n6 criteria: n>=20, pct>=5%, delta>0, FDR<0.1, prop+>=60%, pct_IE>pct_OT",
               fontsize=7.5, fontweight="bold", loc="left")
sz_l = [Line2D([0],[0],marker='o',c='w',markerfacecolor='gray',markersize=np.sqrt(s)/2,alpha=0.4) for s in [25,55,85]]
ax_a.legend(sz_l, ["10%","30%","50%"], fontsize=5.5, loc="upper right", title="ILC-TLS rate", title_fontsize=6, borderpad=0.3)

# Panel B: Dot plot
df_b = df_pos.nlargest(15, "median_delta")
for gi, (_, r) in enumerate(df_b.iterrows()):
    vals = [r["pct_other_TLS"], r["pct_ILC_TLS"]]
    sizes = [15 + max(v, 0.01) * 180 for v in vals]
    colors_dot = ["#4c72b0", "#c44e52"]
    for gj in range(2):
        ax_b.scatter(gj, gi, s=sizes[gj], c=colors_dot[gj], alpha=0.8, edgecolors="black", linewidths=0.3, zorder=5)
ax_b.set_yticks(range(len(df_b)))
ax_b.set_yticklabels(df_b["gene"].values, fontsize=6.5, fontstyle="italic")
ax_b.set_xticks([0, 1])
ax_b.set_xticklabels(["other TLS", "ILC-TLS"], fontsize=6)
ax_b.set_xlim(-0.5, 1.5)
ax_b.set_title(f"b  Detection rate (top {len(df_b)} positive)", fontsize=7.5, fontweight="bold", loc="left")

# Panel C: Paired per-sample for top 4
top4 = [g for g in df_pos.head(4)["gene"].values if g in gene_sample and len(gene_sample[g]) > 5][:4]
for gi, g in enumerate(top4):
    svals = gene_sample[g]
    ie_vals = np.array([v["pct_ie"] for v in svals])
    ot_vals = np.array([v["pct_ot"] for v in svals])
    for i in range(min(len(ie_vals), 100)):
        ax_c.plot([gi-0.25, gi+0.25], [ot_vals[i], ie_vals[i]], color="#ccc", linewidth=0.2, alpha=0.3, zorder=2)
    ax_c.scatter(np.full(len(ot_vals), gi-0.25), ot_vals, s=6, c="#4c72b0", alpha=0.3, zorder=3)
    ax_c.scatter(np.full(len(ie_vals), gi+0.25), ie_vals, s=6, c="#c44e52", alpha=0.3, zorder=3)
    ax_c.scatter(gi-0.25, np.median(ot_vals), s=40, c="#4c72b0", edgecolors="black", linewidths=0.5, zorder=6)
    ax_c.scatter(gi+0.25, np.median(ie_vals), s=40, c="#c44e52", edgecolors="black", linewidths=0.5, zorder=6)
ax_c.set_xticks(range(len(top4)))
ax_c.set_xticklabels(top4, fontsize=7, fontstyle="italic")
ax_c.set_ylabel("pct_expr", fontsize=7)
ax_c.set_title("c  Paired per-sample: other TLS vs ILC-enriched TLS", fontsize=7.5, fontweight="bold", loc="left")
ax_c.legend([Line2D([0],[0],marker='o',c='w',markerfacecolor='#4c72b0',markersize=8),
             Line2D([0],[0],marker='o',c='w',markerfacecolor='#c44e52',markersize=8)],
            ["other TLS", "ILC-TLS"], fontsize=6, loc="upper left")

# Panel D: Category summary
cat_counts = df_pos["category"].value_counts()
cat_colors_plot = [CAT.get(c, "#888") for c in cat_counts.index]
ax_d.barh(range(len(cat_counts)), cat_counts.values, color=cat_colors_plot, alpha=0.6, edgecolor="black", linewidth=0.3)
for i, (cat, n_pos) in enumerate(cat_counts.items()):
    n_total = sum(1 for _, r in df_all.iterrows() if r["category"] == cat)
    ax_d.text(n_pos + 0.3, i, f"{n_pos}/{n_total}", fontsize=6, va="center", color="#666")
ax_d.set_yticks(range(len(cat_counts)))
ax_d.set_yticklabels([c.split("(")[0].strip() for c in cat_counts.index], fontsize=6)
ax_d.set_xlabel("positive receptors", fontsize=7)
ax_d.set_title(f"d  Positive by category ({len(df_pos)} total)", fontsize=7.5, fontweight="bold", loc="left")

fig.tight_layout(pad=0.5)
save_pub(fig, OUT / "fig_niche_positive_4panel")
plt.close()
print(f"\nFigure: {OUT}/fig_niche_positive_4panel.{{svg,pdf,tiff}}")
