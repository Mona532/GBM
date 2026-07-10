"""CellCharter niche-stratified receptor enrichment — all GBM samples"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import lil_matrix
import cellcharter as cc, squidpy as sq, scvi
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

# ── Config ──
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

# ── Niche-stratified per-gene log2FC ──
niche_deltas = []
sample_niches = {}  # for demo figure

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
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct_names.index(c) for c in ILC_TYPES}

        # CellCharter clustering on c2l abundance (all spots)
        for i, c in enumerate(ct_names): adata.obs[c] = q05[:, i]
        try:
            scvi.model.SCVI.setup_anndata(adata, batch_key=None)
            model = scvi.model.SCVI(adata, n_layers=1, n_latent=min(8, len(ct_names)-1))
            model.train(early_stopping=True, max_epochs=20, train_size=0.9)
            adata.obsm['X_scVI'] = model.get_latent_representation()
            sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
            cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter')
            autok = cc.tl.ClusterAutoK(n_clusters=(2, 10), max_runs=3)
            autok.fit(adata, use_rep='X_cellcharter')
            adata.obs['niche'] = autok.predict(adata, use_rep='X_cellcharter')
        except Exception as e:
            print(f"  CellCharter failed for {sd.name}: {e}")
            continue

        niches = sorted(adata.obs['niche'].unique())
        if len(niches) < 2: continue

        # Gene expression (log-normalized)
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr_raw = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values
        lib_size = expr_raw.sum(axis=1) + 1
        expr_norm = np.log1p(expr_raw / (lib_size[:,None]/10000))

        coords_all = adata.obsm["spatial"]

        # For each niche, define ILC-enriched TLS and compare
        sample_niches[sd.name] = {"K": len(niches), "niches": {}}
        for niche_id in niches:
            niche_mask = (adata.obs['niche'] == niche_id).values
            tls_in_niche = tls_mask & niche_mask
            if tls_in_niche.sum() < 6: continue  # need enough TLS in this niche

            # TLS regions within this niche
            tls_gidx = np.where(tls_in_niche)[0]
            tls_coords = coords_all[tls_gidx]
            k_init = min(7, len(tls_gidx))
            knn_all = kneighbors_graph(tls_coords, n_neighbors=k_init, mode="connectivity", include_self=True)
            n_regions, region_labels = connected_components(knn_all, directed=False)

            # Build intra-TLS-region connectivity within this niche
            conn = lil_matrix((len(tls), len(tls)))
            for rid in range(n_regions):
                rm = region_labels == rid
                rg = tls_gidx[rm]
                rc = tls_coords[rm]
                if len(rg) < 3:
                    for gi in rg: conn[gi, gi] = 1
                    continue
                k = min(7, len(rg))
                knn_r = kneighbors_graph(rc, n_neighbors=k, mode="connectivity", include_self=True)
                for li, gi in enumerate(rg):
                    conn[gi, rg[knn_r[li].indices]] = 1
            conn = conn.tocsr()

            # Neighborhood ILC score
            ilc_mat = q05[:, [ilc_idx[c] for c in ILC_TYPES]]
            nhood_ilc = conn @ ilc_mat
            row_sums = np.array(conn.sum(axis=1)).flatten() + 1e-8
            nhood_ilc = nhood_ilc / row_sums[:, None]
            nhood_score = nhood_ilc.max(axis=1)

            p75_nhood = np.percentile(nhood_score[tls_in_niche], 75)
            ilc_enriched = tls_in_niche & (nhood_score >= p75_nhood)
            other_tls = tls_in_niche & ~ilc_enriched

            if ilc_enriched.sum() < 3 or other_tls.sum() < 3: continue

            sample_niches[sd.name]["niches"][niche_id] = {
                "n_tls": int(tls_in_niche.sum()), "n_ilc_enr": int(ilc_enriched.sum()),
                "n_other": int(other_tls.sum()), "n_regions": n_regions
            }

            # Per-gene log2FC within this niche
            ie = expr_norm[ilc_enriched]
            oe = expr_norm[other_tls]
            for g in ALL_RECEPTORS:
                if g not in vn: continue
                gidx = list(vn).index(g)
                im = ie[:, gidx].mean()
                om = oe[:, gidx].mean()
                if im > 0 and om > 0:
                    niche_deltas.append({"sample": sd.name, "niche": niche_id, "gene": g,
                        "log2FC": np.log2(im/om), "category": gene_cat.get(g,"?")})

        print(f"  {sd.name}: {len(niches)} niches, {sum(len(v['niches']) for v in [sample_niches[sd.name]])} with enough TLS")

print(f"\nTotal sample-niche-gene pairs: {len(niche_deltas)}")

# ── Aggregate across samples and niches ──
df_d = pd.DataFrame(niche_deltas)
df_d.to_csv(OUT / "receptor_niche_stratified_raw.csv", index=False)

agg_rows = []
for g in sorted(df_d["gene"].unique()):
    sub = df_d[df_d["gene"] == g]
    if len(sub) < 5: continue
    l2fc = sub["log2FC"].dropna()
    if len(l2fc) < 5: continue
    _, pval = wilcoxon(l2fc, alternative="two-sided")
    agg_rows.append({"gene": g, "n_obs": len(l2fc), "median_log2FC": l2fc.median(),
                     "pvalue": pval, "category": gene_cat.get(g,"?")})

df_agg = pd.DataFrame(agg_rows)
df_agg["fdr"] = false_discovery_control(df_agg["pvalue"].values)
df_agg = df_agg.sort_values("median_log2FC")
df_agg.to_csv(OUT / "receptor_niche_stratified_aggregated.csv", index=False)

# Enriched
df_pos = df_agg[(df_agg["median_log2FC"] > 0.3) & (df_agg["fdr"] < 0.1) & (df_agg["n_obs"] >= 10)]
df_top = df_pos.nlargest(15, "median_log2FC").sort_values("median_log2FC")
print(f"Enriched: {len(df_pos)}, showing top {len(df_top)}")

# ── Figure: niche-stratified receptor enrichment ──
N = len(df_top)
fig, ax = plt.subplots(figsize=(5, N*0.35+0.8))
y = np.arange(N)
for i, (_, r) in enumerate(df_top.iterrows()):
    c = CAT.get(r["category"],"#888")
    ax.plot([0, r["median_log2FC"]], [i, i], color=c, linewidth=1.5, alpha=0.6, solid_capstyle="round", zorder=3)
    ax.scatter(r["median_log2FC"], i, s=55, c=c, alpha=0.9, edgecolors="black", linewidths=0.5, zorder=5)
ax.set_yticks(range(N))
ax.set_yticklabels(df_top["gene"].values, fontsize=8, fontstyle="italic")
ax.set_xlabel("log2(ILC-enriched TLS / other TLS)", fontsize=8, labelpad=6)
ax.set_xlim(-0.2, df_top["median_log2FC"].max()*1.2)
ax.invert_yaxis()
ax.set_title("Receptor enrichment: niche-stratified comparison\n(CellCharter niches, Squidpy intra-TLS neighborhood ILC definition)",
             fontsize=9, fontweight="bold", loc="left", pad=10)
cats_in = set(df_top["category"].values)
leg = [Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in cat_order if c in cats_in]
ax.legend(handles=leg, fontsize=6, loc="lower right", title="Category", title_fontsize=6.5, borderpad=0.3)
n_samples = len(set(d["sample"] for d in niche_deltas))
ax.text(1.0, -0.15, f"n={n_samples} samples, {len(niche_deltas)} niche-gene pairs | CellCharter + Squidpy | Wilcoxon + FDR",
        transform=ax.transAxes, fontsize=5.5, color="#666", ha="right")
fig.tight_layout(pad=1.0)
save_pub(fig, OUT / "fig_receptor_cellcharter_niche")
plt.close()
print(f"Figure: {OUT}/fig_receptor_cellcharter_niche.{{svg,pdf,tiff}}")
print(f"Samples processed: {n_samples}")
