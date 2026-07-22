"""
CellCharter spatial niche → ILC-enriched TLS → receptor comparison (complete, correct)
1. Build AnnData: X = c2l q05 abundance per spot
2. Squidpy spatial_neighbors per sample
3. CellCharter aggregate_neighbors → X_cellcharter (spatially smoothed)
4. CellCharter Cluster → global niche labels
5. Per sample × niche: ILC-enriched TLS vs other TLS
6. Sample-level median delta → Wilcoxon + FDR
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import gc, warnings
warnings.filterwarnings("ignore")

SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): gene_cat[g] = cat_order[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]
CAT = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
OUT = Path(r"E:/GBM/results")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

# ====== Step 1: Build merged AnnData (X = c2l q05, one row per spot) ======
adata_list = []
sample_tls_masks = []  # list of (sample_id, tls_mask)
all_tls_ilc = {c: [] for c in ILC_TYPES}
expr_registry = []  # (sample_id, expr_norm, var_names)

for ds_idx, (tls_dir, h5_dir) in enumerate(DATASETS):
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv).set_index("barcode") if "barcode" in pd.read_csv(tls_csv, nrows=1).columns else pd.read_csv(tls_csv)
        adata = ad.read_h5ad(h5)
        # Align barcodes
        if "barcode" in tls.columns:
            tls = tls.set_index("barcode")
        shared = adata.obs_names.intersection(tls.index)
        if len(shared) < 10: continue
        adata = adata[shared]
        tls = tls.loc[shared]
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 10: continue

        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct_names.index(c) for c in ILC_TYPES}

        # Global TLS ILC
        for c in ILC_TYPES:
            all_tls_ilc[c].extend(q05[tls_mask, ilc_idx[c]].tolist())

        # Gene expression — lib_size from ALL genes, extract only receptor genes
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        rx_present = [g for g in ALL_RECEPTORS if g in vn_all]
        rx_indices = {g: list(vn_all).index(g) for g in rx_present}
        ge_rx = ge[:, [i for g, i in rx_indices.items()]]
        rx_dense = ge_rx.X.toarray() if hasattr(ge_rx.X, "toarray") else np.asarray(ge_rx.X)
        expr_norm = np.log1p(rx_dense / (lib_size[:,None]/10000))
        vn_rx = list(rx_indices.keys())

        # Build mini AnnData for CellCharter (raw q05 abundance preserves cross-sample differences)
        obs_df = pd.DataFrame(index=adata.obs_names)
        obs_df["sample"] = sd.name
        a = ad.AnnData(X=q05.astype(np.float32), obs=obs_df, obsm={"spatial": adata.obsm["spatial"]})
        a.var_names = ct_names

        adata_list.append(a)
        sample_tls_masks.append((sd.name, tls_mask))
        expr_registry.append((sd.name, expr_norm, vn_rx, ilc_idx, q05, rx_indices))

GLOBAL_ILC_THRESH = {c: max(np.percentile(all_tls_ilc[c], 75), 1.0) for c in ILC_TYPES}
print(f"Global ILC thresholds: { {c: round(v,3) for c,v in GLOBAL_ILC_THRESH.items()} }")
print(f"Samples: {len(adata_list)}")

# ====== Step 2: Merge & CellCharter spatial niche ======
adata_all = ad.concat(adata_list, join="inner", label="sample_id")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adata_list])

# Squidpy spatial neighbors (per sample)
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)

# CellCharter: aggregate neighbors → spatial smoothing
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")

# CellCharter clustering
autok = cc.tl.ClusterAutoK(n_clusters=(5, 20), max_runs=5)
autok.fit(adata_all, use_rep="X_cellcharter")
adata_all.obs["niche"] = autok.predict(adata_all, use_rep="X_cellcharter")

best_k = adata_all.obs["niche"].nunique()
print(f"CellCharter niches: {best_k}")
gc.collect()

# ====== Step 3: Per-sample × niche comparison ======
sample_gene_deltas = {}  # gene → [sample_median_delta, ...]

for i, (sid, tls_mask) in enumerate(sample_tls_masks):
    # Get niche labels for this sample
    sample_mask = adata_all.obs["sample"].values == sid
    niche_labels = adata_all.obs["niche"].values[sample_mask]

    # Get expression for this sample
    expr_norm, vn_rx, ilc_idx, q05, rx_idx = expr_registry[i][1:]
    niches_present = sorted(set(niche_labels))
    per_niche_deltas = {}  # gene → [delta_across_niches]

    for niche_id in niches_present:
        niche_mask = niche_labels == niche_id
        tls_in_niche = tls_mask & niche_mask
        n_tls = tls_in_niche.sum()
        if n_tls < 10: continue

        # ILC-enriched TLS: global threshold
        ilc_enr = np.zeros(n_tls, dtype=bool)
        tls_idx = np.where(tls_in_niche)[0]
        for c in ILC_TYPES:
            ilc_enr |= (q05[tls_idx, ilc_idx[c]] >= GLOBAL_ILC_THRESH[c])
        other_tls = ~ilc_enr

        if ilc_enr.sum() < 3 or other_tls.sum() < 3: continue

        ie = expr_norm[tls_idx[ilc_enr]]
        oe = expr_norm[tls_idx[other_tls]]
        for gi, g in enumerate(vn_rx):
            delta = ie[:, gi].mean() - oe[:, gi].mean()
            per_niche_deltas.setdefault(g, []).append(delta)

    # Sample-level: median delta across niches
    for g, deltas in per_niche_deltas.items():
        if len(deltas) == 0: continue
        sample_gene_deltas.setdefault(g, []).append(np.median(deltas))

    if (i+1) % 20 == 0: print(f"  Processed {i+1}/{len(sample_tls_masks)} samples")

# Free expression matrices from memory
expr_registry.clear()
gc.collect()

# ====== Step 4: Cross-sample statistics ======
rows = []
for g, sample_deltas in sample_gene_deltas.items():
    arr = np.array(sample_deltas)
    if len(arr) < 5: continue
    try:
        _, pval = wilcoxon(arr, alternative="two-sided")
    except:
        pval = np.nan
    if np.isnan(pval) or pval <= 0 or pval >= 1: continue
    rows.append({"gene": g, "n_samples": len(arr), "median_delta": np.median(arr),
                 "pvalue": pval, "category": gene_cat.get(g,"?")})

df = pd.DataFrame(rows)
if len(df) == 0:
    print("No genes passed the filtering criteria. Check thresholds.")
    exit(0)
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("median_delta")
df.to_csv(OUT / "receptor_cellcharter_final.csv", index=False)

df_pos = df[(df["median_delta"] > 0.0005) & (df["fdr"] < 0.1)]
df_top = df_pos.nlargest(15, "median_delta").sort_values("median_delta")
n_top = len(df_top)
print(f"\nEnriched genes FDR<0.1: {len(df_pos)}, showing top {n_top}")

for _, r in df_top.iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} Δ={r['median_delta']:+.5f} FDR={r['fdr']:.4f} n={r['n_samples']} {s}")

# ====== Figure ======
if n_top > 0:
    N = n_top
    fig, ax = plt.subplots(figsize=(5, N*0.38+1))
    y = np.arange(N)
    for i, (_, r) in enumerate(df_top.iterrows()):
        c = CAT.get(r["category"],"#888")
        ax.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.6, solid_capstyle="round", zorder=3)
        ax.scatter(r["median_delta"], i, s=55, c=c, alpha=0.9, edgecolors="black", linewidths=0.5, zorder=5)
    ax.set_yticks(range(N))
    ax.set_yticklabels(df_top["gene"].values, fontsize=8, fontstyle="italic")
    ax.set_xlabel("Δ mean log-expr (ILC-enriched TLS − other TLS)", fontsize=8, labelpad=6)
    ax.set_xlim(-0.01, df_top["median_delta"].max()*1.2)
    ax.invert_yaxis()
    ax.set_title("Receptor expression: niche-stratified\n(CellCharter on c2l abundance, Squidpy spatial graph)", fontsize=9, fontweight="bold", loc="left", pad=10)
    cats_in = set(df_top["category"].values)
    ax.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in cat_order if c in cats_in],
              fontsize=6, loc="lower right", title="Category", title_fontsize=6.5, borderpad=0.3)
    ax.text(1.0, -0.15, f"n={df_top['n_samples'].max()} samples | CellCharter | Squidpy | niche-stratified | Wilcoxon+FDR",
            transform=ax.transAxes, fontsize=5.5, color="#666", ha="right")
    fig.tight_layout(pad=1.0)
    save_pub(fig, OUT / "fig_receptor_cellcharter_final")
    plt.close()
    print(f"Figure: {OUT}/fig_receptor_cellcharter_final.{{svg,pdf,tiff}}")
else:
    print("No enriched genes to plot.")
