"""
CellCharter niche-stratified receptor analysis — all 6 fixes applied
Input:  c2l abundance (23 cell types), NOT gene expression
Niche:  joint PCA + GMM across all 142 samples → global niche labels
ILC:    global TLS P75 threshold, >= 1.0, n_tls_in_niche >= 20
Graph:  Squidpy spatial_neighbors (distance-based)
Effect: delta = mean_log_expr_diff (ILC-enriched − other TLS)
Stats:  sample-level median delta → Wilcoxon + FDR across samples
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

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

# ====== Step 1: Collect c2l abundance, gene expression, global TLS ILC ======
all_q05_scaled = []
sample_records = []  # (sid, n_spots, q05_scaled, q05_raw, ilc_idx, tls_mask, expr_norm, var_names, coords)
global_tls_ilc = {c: [] for c in ILC_TYPES}

for tls_dir, h5_dir in DATASETS:
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
        # Gene expression
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr_raw = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values
        lib_size = expr_raw.sum(axis=1) + 1
        expr_norm = np.log1p(expr_raw / (lib_size[:,None]/10000))
        # Global TLS ILC
        for c in ILC_TYPES:
            global_tls_ilc[c].extend(q05[tls_mask, ilc_idx[c]].tolist())
        # Standardize c2l per sample
        scaler = StandardScaler()
        q05_scaled = scaler.fit_transform(q05)
        all_q05_scaled.append(q05_scaled)
        coords = adata.obsm["spatial"]
        sample_records.append((sd.name, q05.shape[0], q05_scaled, q05, ilc_idx, tls_mask, expr_norm, vn, coords))

GLOBAL_ILC_THRESH = {c: max(np.percentile(global_tls_ilc[c], 75), 1.0) for c in ILC_TYPES}
print(f"Global ILC thresholds: { {c: round(v,3) for c,v in GLOBAL_ILC_THRESH.items()} }")
print(f"Samples: {len(sample_records)}")

# ====== Step 2: Joint PCA + GMM → global niche labels ======
stacked = np.vstack(all_q05_scaled)
pca = PCA(n_components=15)
pca_data = pca.fit_transform(stacked)
print(f"PCA: 15 PCs, explained var={pca.explained_variance_ratio_.sum():.2f}")

n_sample = min(100000, pca_data.shape[0])
idx = np.random.RandomState(42).choice(pca_data.shape[0], n_sample, replace=False)
bics = [GaussianMixture(n_components=k, covariance_type='tied', random_state=42).fit(pca_data[idx]).bic(pca_data[idx]) for k in range(2, 16)]
best_k = np.argmin(bics) + 2
print(f"Best K = {best_k}")

gmm = GaussianMixture(n_components=best_k, covariance_type='tied', random_state=42)
gmm.fit(pca_data[idx])
offset = 0
niche_data = []  # (sid, niche_labels, q05, ilc_idx, tls_mask, expr_norm, vn, coords)
for sid, n_spots, qs, qr, ilc_idx, tls_mask, expr_norm, vn, coords in sample_records:
    labels = gmm.predict(pca_data[offset:offset+n_spots])
    niche_data.append((sid, labels, qr, ilc_idx, tls_mask, expr_norm, vn, coords))
    offset += n_spots

print(f"Niches: {best_k}, samples assigned")

# ====== Step 3: Per-niche delta (ILC-enriched TLS vs other TLS) ======
sample_gene_deltas = {}  # gene -> list of (sample, median_delta)

for sid, niche_labels, q05, ilc_idx, tls_mask, expr_norm, vn, coords in niche_data:
    niches_present = sorted(set(niche_labels))
    niche_gene_deltas = {}  # gene -> list of deltas across niches in this sample

    for niche_id in niches_present:
        niche_mask = niche_labels == niche_id
        tls_in_niche = tls_mask & niche_mask
        n_tls = tls_in_niche.sum()
        if n_tls < 20: continue

        # Squidpy spatial neighbors on niche spots, then restrict to TLS
        temp_adata = ad.AnnData(X=np.zeros((len(tls_mask), 1)), obsm={"spatial": coords})
        try:
            sq.gr.spatial_neighbors(temp_adata, coord_type="generic", delaunay=True)
            conn_full = temp_adata.obsp["spatial_connectivities"]
        except:
            from sklearn.neighbors import kneighbors_graph
            conn_full = kneighbors_graph(coords, n_neighbors=7, mode="connectivity", include_self=True)

        # Neighborhood ILC: only using connections within TLS-in-niche
        tls_idx = np.where(tls_in_niche)[0]
        if len(tls_idx) < 3: continue
        # Build subgraph: only edges where BOTH ends are TLS-in-niche
        conn_tls = conn_full[tls_idx][:, tls_idx]
        ilc_mat = q05[:, [ilc_idx[c] for c in ILC_TYPES]]
        ilc_mat_tls = ilc_mat[tls_idx]
        nhood_ilc = conn_tls @ ilc_mat_tls
        row_sums = np.array(conn_tls.sum(axis=1)).flatten() + 1e-8
        nhood_ilc = nhood_ilc / row_sums[:, None]
        nhood_score = nhood_ilc.max(axis=1)

        # ILC-enriched: global threshold (not per-niche P75)
        ilc_enr_local = np.zeros(len(tls_idx), dtype=bool)
        for c in ILC_TYPES:
            ilc_enr_local |= (q05[tls_idx, ilc_idx[c]] >= GLOBAL_ILC_THRESH[c])
        # Intersect with top ILC nhood score
        p75_nhood = np.percentile(nhood_score, 75)
        ilc_enriched = ilc_enr_local & (nhood_score >= p75_nhood)
        other_tls = np.ones(len(tls_idx), dtype=bool) & ~ilc_enriched

        if ilc_enriched.sum() < 3 or other_tls.sum() < 3: continue

        # Per-gene mean diff (delta)
        ie = expr_norm[tls_idx[ilc_enriched]]
        oe = expr_norm[tls_idx[other_tls]]
        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)
            delta = ie[:, gidx].mean() - oe[:, gidx].mean()
            niche_gene_deltas.setdefault(g, []).append(delta)

    # Sample-level aggregate: median delta per gene
    for g, deltas in niche_gene_deltas.items():
        if len(deltas) == 0: continue
        sample_gene_deltas.setdefault(g, []).append(np.median(deltas))

# ====== Step 4: Cross-sample Wilcoxon test ======
rows = []
for g, sample_deltas in sample_gene_deltas.items():
    arr = np.array(sample_deltas)
    if len(arr) < 5: continue
    _, pval = wilcoxon(arr, alternative="two-sided")
    rows.append({"gene": g, "n_samples": len(arr), "median_delta": np.median(arr),
                 "pvalue": pval, "category": gene_cat.get(g,"?")})

df = pd.DataFrame(rows)
df = df[df["pvalue"].notna() & (df["pvalue"] > 0) & (df["pvalue"] < 1)]
df["fdr"] = false_discovery_control(df["pvalue"].values)
df = df.sort_values("median_delta")
df.to_csv(OUT / "receptor_niche_stratified_v2.csv", index=False)

df_pos = df[df["median_delta"] > 0.0005][df["fdr"] < 0.1]
df_top = df_pos.nlargest(15, "median_delta").sort_values("median_delta")
print(f"\nEnriched: {len(df_pos)}, showing top {len(df_top)}")

for _, r in df_top.iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s} Δ={r['median_delta']:+.5f} FDR={r['fdr']:.4f} n={r['n_samples']} {s}")

# ====== Figure ======
N = len(df_top)
if N > 0:
    fig, ax = plt.subplots(figsize=(5, N*0.38+1))
    y = np.arange(N)
    for i, (_, r) in enumerate(df_top.iterrows()):
        c = CAT.get(r["category"],"#888")
        ax.plot([0, r["median_delta"]], [i, i], color=c, linewidth=1.5, alpha=0.6, solid_capstyle="round", zorder=3)
        ax.scatter(r["median_delta"], i, s=55, c=c, alpha=0.9, edgecolors="black", linewidths=0.5, zorder=5)
    ax.set_yticks(range(N))
    ax.set_yticklabels(df_top["gene"].values, fontsize=8, fontstyle="italic")
    ax.set_xlabel("Δ mean log-expr (ILC-enriched − other TLS)", fontsize=8, labelpad=6)
    ax.set_xlim(-0.01, df_top["median_delta"].max()*1.2)
    ax.invert_yaxis()
    ax.set_title("Receptor expression: niche-stratified comparison\n(CellCharter on c2l abundance, global niches)", fontsize=9, fontweight="bold", loc="left", pad=10)
    cats_in = set(df_top["category"].values)
    ax.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in cat_order if c in cats_in],
              fontsize=6, loc="lower right", title="Category", title_fontsize=6.5, borderpad=0.3)
    ax.text(1.0, -0.15, f"n={df_top['n_samples'].max()} samples | CellCharter c2l niches | Squidpy spatial graph | Wilcoxon + FDR",
            transform=ax.transAxes, fontsize=5.5, color="#666", ha="right")
    fig.tight_layout(pad=1.0)
    save_pub(fig, OUT / "fig_receptor_cellcharter_v2")
    plt.close()
    print(f"Figure: {OUT}/fig_receptor_cellcharter_v2.{{svg,pdf,tiff}}")
else:
    print("No enriched genes found.")
