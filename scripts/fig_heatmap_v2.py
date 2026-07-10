"""
Heatmap v2: Informative niche composition + ILC-rich TLS distribution + receptor/ligand detection
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import wilcoxon, false_discovery_control
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})
def save_pub(fig, stem):
    for fmt in ["svg","pdf","tiff"]: fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=600)

CAT = {"Excit (Glutamate)":"#E64A19","Inhib (GABA/Gly)":"#2E7D32","Cholinergic (ACh)":"#1565C0","DA/NE":"#7B1FA2","Serotonin (5-HT)":"#C62828"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
RECEPTOR_GENES = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): RECEPTOR_GENES[g] = list(CAT.keys())[idx]
LIGANDS = ["NMU","VIP","ADM","CALCA","CALCB","NPY","TAC1","TAC3","NTS","CCK","GRP","GAL","PENK","PDYN","PNOC","CRH","UCN","AGRP","POMC","MCH"]
ALL_GENES = sorted(set(list(RECEPTOR_GENES.keys()) + LIGANDS))
ILC_TYPES = ["ILC1","ILC2","ILC3"]
GLOBAL_THRESH = {"ILC1":1.034,"ILC2":1.0,"ILC3":1.035}
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results")
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# ====== Step 1: Global CellCharter ======
adata_list, sample_data = [], []
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
        ct_names = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx_map = {c: ct_names.index(c) for c in ILC_TYPES}
        coords = adata.obsm["spatial"]
        q05_z = (q05 - q05.mean(axis=0)) / (q05.std(axis=0) + 1e-8)
        obs_df = pd.DataFrame(index=adata.obs_names); obs_df["sample"] = sd.name
        a = ad.AnnData(X=q05_z.astype(np.float32), obs=obs_df, obsm={"spatial": coords})
        a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
        adata_list.append(a)
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        vn_all = ge.var_names.values
        genes_present = [g for g in ALL_GENES if g in vn_all]
        gene_idx = [list(vn_all).index(g) for g in genes_present]
        rx_raw = ge[:, gene_idx].X.toarray() if hasattr(ge[:, gene_idx].X,"toarray") else np.asarray(ge[:, gene_idx].X)
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))
        sample_data.append((sd.name, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords))

adata_all = ad.concat(adata_list, join="inner")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adata_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
# Fixed K=5 for reproducibility
gmm = cc.tl.Cluster(n_clusters=5, random_state=42)
gmm.fit(adata_all, use_rep="X_cellcharter")
niche_all = gmm.predict(adata_all, use_rep="X_cellcharter")
K = len(set(niche_all))
print(f"Global niches: {K}")

# ====== Step 2: Niche composition (log2 enrichment over global mean) ======
niche_all_int = niche_all.astype(int)
# Global mean per cell type
all_q05_z = np.vstack([sd[3] for sd in sample_data])
global_mean = all_q05_z.mean(axis=0)  # ~0 since z-scored, use raw instead

# Use raw q05 for meaningful enrichment
all_q05_raw = np.vstack([sd[2] for sd in sample_data])
global_mean_raw = all_q05_raw.mean(axis=0) + 1e-8
niche_composition = np.zeros((K, len(ct_names)))
for ni in range(K):
    mask = niche_all_int == ni
    if mask.sum() > 0:
        niche_mean = all_q05_raw[mask].mean(axis=0)
        # log2 enrichment over global
        niche_composition[ni] = np.log2(niche_mean / global_mean_raw)

# ILC-rich TLS data
gene_sample = {}
niche_ilc_counts_all = {ni: 0 for ni in range(K)}
offset = 0
for sid, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords in sample_data:
    n = len(tls_mask); niches = niche_all_int[offset:offset+n]; offset += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
    ilc_rich = ilc_high & tls_mask; other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue
    for ni in range(K):
        niche_ilc_counts_all[ni] += int((niches == ni).sum())
    ie = expr_norm[ilc_rich]; oe = expr_norm[other_tls]
    for gi, g in enumerate(genes_present):
        pct_ilc = (rx_raw[ilc_rich, gi] > 0).mean()
        pct_other = (rx_raw[other_tls, gi] > 0).mean()
        gene_sample.setdefault(g, []).append({"pct_ilc":pct_ilc,"pct_other":pct_other,"delta":ie[:,gi].mean()-oe[:,gi].mean()})

# ILC-rich enrichment per niche: % of ILC-rich TLS spots in each niche
ilc_total = sum(niche_ilc_counts_all.values()) + 1e-8
ilc_niche_pct = {ni: niche_ilc_counts_all[ni]/ilc_total*100 for ni in range(K)}

# Aggregate receptor/ligand
rows = []
for g, svals in gene_sample.items():
    n_s = len(svals); pct_i = np.array([v["pct_ilc"] for v in svals])
    pct_o = np.array([v["pct_other"] for v in svals]); deltas = np.array([v["delta"] for v in svals])
    rows.append({"gene":g,"n":n_s,"pct_ILC":np.median(pct_i),"pct_other":np.median(pct_o),
        "median_delta":np.median(deltas),"type":"receptor" if g in RECEPTOR_GENES else "ligand",
        "category":RECEPTOR_GENES.get(g,"Ligand")})
df_all = pd.DataFrame(rows).sort_values("pct_ILC", ascending=False)
df_all.to_csv(OUT / "receptor_ligand_simple.csv", index=False)

# ====== Figure ======
fig = plt.figure(figsize=(13, 9))

# --- Panel A: Niche composition (log2 enrichment) ---
ax_a = fig.add_axes([0.04, 0.53, 0.44, 0.42])
vmax = max(abs(niche_composition.min()), abs(niche_composition.max()), 1.0)
im_a = ax_a.imshow(niche_composition, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax_a.set_xticks(range(len(ct_names)))
ax_a.set_xticklabels(ct_names, fontsize=4.5, rotation=90)
ax_a.set_yticks(range(K))
ax_a.set_yticklabels([f"N{n}" for n in range(K)], fontsize=5.5)
# ILC-rich % annotation on right
for ni in range(K):
    pct = ilc_niche_pct.get(ni, 0)
    ax_a.text(len(ct_names), ni, f" {pct:.1f}%", fontsize=5, color="#c44e52" if pct > 10 else "#666", va="center")
ax_a.set_title(f"a  Niche composition (log2 enrichment over global mean)\nRight: % ILC-rich TLS spots, {K} niches", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_a, ax=ax_a, shrink=0.6)

# --- Panel B: ILC-rich TLS enrichment per niche ---
ax_b = fig.add_axes([0.54, 0.53, 0.43, 0.42])
niches_sorted = sorted(ilc_niche_pct.keys(), key=lambda n: ilc_niche_pct[n], reverse=True)
pcts = [ilc_niche_pct[n] for n in niches_sorted]
colors_b = ["#c44e52" if p > 10 else "#4c72b0" if p > 5 else "#aaa" for p in pcts]
bars = ax_b.bar(range(K), pcts, color=colors_b, alpha=0.7, edgecolor="black", linewidth=0.3)
ax_b.set_xticks(range(K))
ax_b.set_xticklabels([f"N{n}" for n in niches_sorted], fontsize=6)
ax_b.set_ylabel("% ILC-rich TLS spots", fontsize=7)
ax_b.set_title(f"b  ILC-rich TLS distribution across niches", fontsize=7.5, fontweight="bold", loc="left")
ax_b.axhline(y=100/K, color="#ccc", linestyle="--", linewidth=0.5, label=f"uniform ({100/K:.0f}%)")
ax_b.legend(fontsize=6)

# --- Panel C: Receptor detection ---
df_rcp = df_all[df_all["type"]=="receptor"].head(30)
rcp_mat = df_rcp[["pct_other","pct_ILC"]].values
ax_c = fig.add_axes([0.04, 0.04, 0.30, 0.40])
im_c = ax_c.imshow(rcp_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_c.set_xticks([0, 1]); ax_c.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_c.set_yticks(range(len(df_rcp))); ax_c.set_yticklabels(df_rcp["gene"].values, fontsize=5.5, fontstyle="italic")
ax_c.set_title(f"c  Receptor detection rate", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_c, ax=ax_c, shrink=0.7)

# --- Panel D: Ligand detection ---
df_lig = df_all[df_all["type"]=="ligand"]
lig_mat = df_lig[["pct_other","pct_ILC"]].values
ax_d = fig.add_axes([0.37, 0.04, 0.20, 0.40])
im_d = ax_d.imshow(lig_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_d.set_xticks([0, 1]); ax_d.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_d.set_yticks(range(len(df_lig))); ax_d.set_yticklabels(df_lig["gene"].values, fontsize=5.5, fontstyle="italic")
ax_d.set_title(f"d  Ligand detection", fontsize=7.5, fontweight="bold", loc="left")
plt.colorbar(im_d, ax=ax_d, shrink=0.7)

# --- Panel E: Category legend ---
ax_e = fig.add_axes([0.63, 0.04, 0.34, 0.40])
ax_e.axis("off")
# Show receptor pct_ILC ranking with category colors
for j, (_, r) in enumerate(df_rcp.head(25).iterrows()):
    c = CAT.get(r["category"], "#888")
    pct = r["pct_ILC"]
    ax_e.barh(j, pct, color=c, alpha=0.6, height=0.7)
    ax_e.text(pct + 0.01, j, r["gene"], fontsize=5.5, fontstyle="italic", va="center")
ax_e.set_ylim(25.5, -0.5)
ax_e.set_xlim(0, 0.6)
ax_e.set_title("e  Receptor pct_expr (ILC-rich TLS)", fontsize=7.5, fontweight="bold", loc="left")
ax_e.set_xlabel("pct_expr", fontsize=6)

fig.savefig(OUT / "fig_heatmap_v2.svg", bbox_inches="tight")
fig.savefig(OUT / "fig_heatmap_v2.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_heatmap_v2.tiff", dpi=600, bbox_inches="tight")
plt.close()
print(f"Figure: {OUT}/fig_heatmap_v2.{{svg,pdf,tiff}}")
print(f"Niche composition range: [{niche_composition.min():.2f}, {niche_composition.max():.2f}]")
for ni in range(K):
    top_ct = ct_names[np.argmax(niche_composition[ni])]
    print(f"  N{ni}: top={top_ct} ({niche_composition[ni].max():.2f}), ILC-rich={ilc_niche_pct.get(ni,0):.1f}%")
