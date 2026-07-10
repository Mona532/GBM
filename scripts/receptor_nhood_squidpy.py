"""Squidpy spatial neighborhood analysis + receptor enrichment with Nature figures"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
from scipy.sparse.csgraph import connected_components
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import squidpy as sq
import warnings
warnings.filterwarnings("ignore")

# ── Nature defaults ──
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.5, "legend.frameon": False,
})

CAT_COLORS = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}

rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): gene_cat[g] = cat_order[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
OUT = Path(r"E:/GBM/results")

def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]:
        fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

# ====== Step 1: Spatial neighbors, neighborhood ILC per spot ======
nhood_deltas = []
demo_sample = None

for ds_idx, (tls_dir, h5_dir) in enumerate(DATASETS):
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 10: continue
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}

        # Save one sample for demo figure
        if demo_sample is None and ds_idx == 0:
            demo_sample = (sd.name, adata, q05, ct, ilc_idx, tls, tls_mask)

        # Identify individual TLS regions (connected components within TLS spots)
        coords_all = adata.obsm["spatial"]
        tls_global_idx = np.where(tls_mask)[0]
        tls_coords = coords_all[tls_mask]
        # Full k-NN on all TLS to find which are connected
        from sklearn.neighbors import kneighbors_graph
        k_conn = min(7, tls_mask.sum())
        knn_all_tls = kneighbors_graph(tls_coords, n_neighbors=k_conn, mode="connectivity", include_self=True)
        # Find connected components = individual TLS regions
        from scipy.sparse.csgraph import connected_components
        n_regions, region_labels = connected_components(knn_all_tls, directed=False)
        print(f"  {sd.name}: {n_regions} TLS regions, sizes: {np.bincount(region_labels)}")

        # Build neighborhood within each TLS region separately
        from scipy.sparse import lil_matrix
        conn = lil_matrix((len(tls), len(tls)))
        for region_id in range(n_regions):
            region_mask = region_labels == region_id
            region_global = tls_global_idx[region_mask]
            region_coords = tls_coords[region_mask]
            if len(region_global) < 3:
                # tiny region: just add self-connections
                for gi in region_global:
                    conn[gi, gi] = 1
                continue
            k = min(7, len(region_global))
            knn_region = kneighbors_graph(region_coords, n_neighbors=k, mode="connectivity", include_self=True)
            for local_i, global_i in enumerate(region_global):
                conn[global_i, region_global[knn_region[local_i].indices]] = 1
        conn = conn.tocsr()

        # Neighborhood ILC within TLS: only average over neighbor TLS spots
        ilc_mat = q05[:, [ilc_idx[c] for c in ILC_TYPES]]
        nhood_ilc = conn @ ilc_mat  # only TLS spots get neighborhood values
        row_sums = np.array(conn.sum(axis=1)).flatten()
        row_sums[row_sums==0] = 1
        nhood_ilc = nhood_ilc / row_sums[:, None]

        # ILC score: max across 3 ILCs
        nhood_ilc_score = nhood_ilc.max(axis=1)

        # Global P75 of neighborhood ILC score within TLS spots
        p75_nhood = np.percentile(nhood_ilc_score[tls_mask], 75)

        # ILC-enriched TLS (neighborhood definition)
        ilc_enriched = tls_mask & (nhood_ilc_score >= p75_nhood)
        other_tls = tls_mask & ~ilc_enriched

        if ilc_enriched.sum() < 3 or other_tls.sum() < 3:
            continue

        # Gene expression
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr_raw = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values
        lib_size = expr_raw.sum(axis=1)
        expr_norm = np.log1p(expr_raw / (lib_size[:,None]/10000 + 1))

        ilc_expr = expr_norm[ilc_enriched]
        oth_expr = expr_norm[other_tls]
        non_expr = expr_norm[~tls_mask]

        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)
            ile = ilc_expr[:, gidx].mean()
            ote = oth_expr[:, gidx].mean()
            if ile > 0 and ote > 0:
                nhood_deltas.append({"sample":sd.name,"gene":g,
                    "log2FC_nhood":np.log2(ile/ote),"mean_ILC_enriched":ile,"mean_other_TLS":ote,
                    "category":gene_cat.get(g,"?")})

print(f"Samples with ILC-enriched TLS: {len(set(d['sample'] for d in nhood_deltas))}")

# ====== Step 2: Demo figure on one sample ======
if demo_sample:
    sname, adata_d, q05_d, ct_d, ilc_idx_d, tls_d, tls_mask_d = demo_sample
    # Build neighbors within TLS regions only
    coords_d = adata_d.obsm["spatial"]
    tls_gidx = np.where(tls_mask_d)[0]
    tls_coords_d = coords_d[tls_gidx]
    from sklearn.neighbors import kneighbors_graph
    k = min(7, len(tls_gidx))
    knn_d = kneighbors_graph(tls_coords_d, n_neighbors=k, mode="connectivity", include_self=True)
    n_regions, region_labels = connected_components(knn_d, directed=False)
    from scipy.sparse import lil_matrix
    conn_d = lil_matrix((len(tls_d), len(tls_d)))
    for rid in range(n_regions):
        rm = region_labels == rid
        rg = tls_gidx[rm]; rc = tls_coords_d[rm]
        if len(rg) < 3:
            for gi in rg: conn_d[gi, gi] = 1
            continue
        k2 = min(7, len(rg))
        knn_r = kneighbors_graph(rc, n_neighbors=k2, mode="connectivity", include_self=True)
        for li, gi in enumerate(rg):
            conn_d[gi, rg[knn_r[li].indices]] = 1
    conn_d = conn_d.tocsr()
    ilc_mat_d = q05_d[:, [ilc_idx_d[c] for c in ILC_TYPES]]
    nhood_ilc_d = conn_d @ ilc_mat_d
    row_sums = np.array(conn_d.sum(axis=1)).flatten() + 1e-8
    nhood_ilc_d = nhood_ilc_d / row_sums[:, None]
    nhood_score = nhood_ilc_d.max(axis=1)
    p75_nh = np.percentile(nhood_score[tls_mask_d], 75)
    ilc_enr = tls_mask_d & (nhood_score >= p75_nh)

    coords = adata_d.obsm["spatial"]
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))

    # Panel A: TLS mask
    ax = axes[0]
    ax.scatter(coords[:,0], coords[:,1], c="lightgray", s=0.5, rasterized=True)
    ax.scatter(coords[tls_mask_d,0], coords[tls_mask_d,1], c="#c44e52", s=2, rasterized=True)
    ax.set_title(f"{sname}\nTLS spots", fontsize=7, fontweight="bold")
    ax.axis("off")

    # Panel B: Neighborhood ILC score
    ax = axes[1]
    s = ax.scatter(coords[:,0], coords[:,1], c=nhood_score, s=0.5, cmap="YlOrRd", rasterized=True)
    ax.scatter(coords[tls_mask_d,0], coords[tls_mask_d,1], s=3, facecolors="none", edgecolors="red", linewidths=0.3)
    ax.set_title(f"Neighborhood ILC score\n(P75={p75_nh:.2f})", fontsize=7, fontweight="bold")
    ax.axis("off")
    plt.colorbar(s, ax=ax, shrink=0.6)

    # Panel C: ILC-enriched TLS
    ax = axes[2]
    ax.scatter(coords[:,0], coords[:,1], c="lightgray", s=0.5, rasterized=True)
    ax.scatter(coords[other_tls := tls_mask_d & ~ilc_enr,0], coords[other_tls,1], c="#4c72b0", s=2, label="Other TLS")
    ax.scatter(coords[ilc_enr,0], coords[ilc_enr,1], c="#c44e52", s=5, edgecolors="black", linewidths=0.3, label="ILC-enriched TLS")
    ax.set_title(f"ILC-enriched TLS (n={ilc_enr.sum()})\nvs other TLS (n={other_tls.sum()})", fontsize=7, fontweight="bold")
    ax.axis("off")
    ax.legend(fontsize=6, loc="lower left")

    fig.tight_layout(pad=0.5)
    save_pub(fig, OUT / "fig_nhood_ilc_demo")
    plt.close()
    print(f"Demo figure: {OUT}/fig_nhood_ilc_demo.{{svg,pdf,tiff}}")

# ====== Step 3: Aggregate receptor enrichment ======
df_d = pd.DataFrame(nhood_deltas)
agg_rows = []
for g in sorted(df_d["gene"].unique()):
    sub = df_d[df_d["gene"]==g]
    if len(sub) < 5: continue
    l2fc = sub["log2FC_nhood"].dropna()
    if len(l2fc) < 5: continue
    _, pval = wilcoxon(l2fc, alternative="two-sided")
    agg_rows.append({"gene":g,"n_samples":len(l2fc),"median_log2FC":l2fc.median(),"pvalue":pval,"category":gene_cat.get(g,"?")})

df_agg = pd.DataFrame(agg_rows)
df_agg["fdr"] = false_discovery_control(df_agg["pvalue"].values)
df_agg = df_agg.sort_values("median_log2FC")
df_agg.to_csv(OUT / "receptor_nhood_enrichment.csv", index=False)

# Enriched genes
df_pos = df_agg[(df_agg["median_log2FC"]>0.3)&(df_agg["fdr"]<0.1)&(df_agg["n_samples"]>=10)]
df_top = df_pos.nlargest(12, "median_log2FC").sort_values("median_log2FC")
N = len(df_top)
print(f"Enriched genes: {len(df_pos)}, showing top {N}")

# Figure
fig, ax = plt.subplots(figsize=(5, N*0.35+0.8))
y = np.arange(N)
for i, (_, r) in enumerate(df_top.iterrows()):
    c = CAT_COLORS.get(r["category"],"#888")
    ax.plot([0, r["median_log2FC"]], [i, i], color=c, linewidth=1.5, alpha=0.6, solid_capstyle="round", zorder=3)
    ax.scatter(r["median_log2FC"], i, s=55, c=c, alpha=0.9, edgecolors="black", linewidths=0.5, zorder=5)
ax.set_yticks(range(N))
ax.set_yticklabels(df_top["gene"].values, fontsize=8, fontstyle="italic")
ax.set_xlabel("log2(ILC-enriched TLS / other TLS)", fontsize=8, labelpad=6)
ax.set_xlim(-0.2, df_top["median_log2FC"].max()*1.2)
ax.invert_yaxis()
ax.set_title("Receptor expression: ILC-enriched vs other TLS\n(Squidpy neighborhood ILC definition)", fontsize=9, fontweight="bold", loc="left", pad=10)

cats_present = set(df_top["category"].values)
from matplotlib.patches import Patch
leg = [Patch(color=CAT_COLORS[c], alpha=0.5, label=c.split("(")[0].strip()) for c in cat_order if c in cats_present]
ax.legend(handles=leg, fontsize=6, loc="lower right", title="Category", title_fontsize=6.5, borderpad=0.3)
ax.text(1.0, -0.15, f"n = {df_top['n_samples'].max()} samples | Squidpy spatial neighbors (k=7) | Wilcoxon + FDR",
        transform=ax.transAxes, fontsize=5.5, color="#666", ha="right")

fig.tight_layout(pad=1.0)
save_pub(fig, OUT / "fig_receptor_nhood_enrichment")
plt.close()
print(f"Figure: {OUT}/fig_receptor_nhood_enrichment.{{svg,pdf,tiff}}")
print(f"CSV: {OUT}/receptor_nhood_enrichment.csv")
