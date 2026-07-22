"""Nature publication figure: niche composition, ILC-rich TLS distribution, receptor/ligand detection"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

# ── Nature defaults ──
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.5,
    "legend.frameon": False, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
})
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

# ====== Data collection ======
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

# ====== CellCharter K=5 ======
adata_all = ad.concat(adata_list, join="inner")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adata_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
gmm = cc.tl.Cluster(n_clusters=5, random_state=42)
gmm.fit(adata_all, use_rep="X_cellcharter")
niche_all = gmm.predict(adata_all, use_rep="X_cellcharter")
K = 5

# ====== Niche composition ======
all_q05_raw = np.vstack([sd[2] for sd in sample_data])
global_mean_raw = all_q05_raw.mean(axis=0) + 1e-8
niche_comp = np.zeros((K, len(ct_names)))
for ni in range(K):
    mask = niche_all == ni
    if mask.sum() > 0:
        niche_comp[ni] = np.log2(all_q05_raw[mask].mean(axis=0) / global_mean_raw)

# ====== ILC-rich TLS data ======
gene_sample = {}; ilc_niche_counts = np.zeros(K)
offset = 0
for sid, tls_mask, q05, q05_z, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords in sample_data:
    n = len(tls_mask); niches = niche_all[offset:offset+n]; offset += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx_map[c]] >= GLOBAL_THRESH[c])
    ilc_rich = ilc_high & tls_mask; other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue
    for ni in range(K): ilc_niche_counts[ni] += int((niches == ni).sum())
    ie = expr_norm[ilc_rich]; oe = expr_norm[other_tls]
    for gi, g in enumerate(genes_present):
        pct_ilc = (rx_raw[ilc_rich, gi] > 0).mean(); pct_other = (rx_raw[other_tls, gi] > 0).mean()
        gene_sample.setdefault(g, []).append({"pct_ilc":pct_ilc,"pct_other":pct_other})

rows = []
for g, svals in gene_sample.items():
    pct_i = np.array([v["pct_ilc"] for v in svals]); pct_o = np.array([v["pct_other"] for v in svals])
    rows.append({"gene":g,"pct_ILC":np.median(pct_i),"pct_other":np.median(pct_o),
        "type":"receptor" if g in RECEPTOR_GENES else "ligand","category":RECEPTOR_GENES.get(g,"Ligand")})
df_all = pd.DataFrame(rows)

ilc_total = ilc_niche_counts.sum() + 1e-8
ilc_niche_pct = {ni: ilc_niche_counts[ni]/ilc_total*100 for ni in range(K)}

# ====== Figure ======
fig = plt.figure(figsize=(10, 11))

# ── Panel A: Niche composition ──
ax_a = fig.add_axes([0.06, 0.68, 0.42, 0.28])
vmax = max(abs(niche_comp.min()), abs(niche_comp.max()), 1.0)
im_a = ax_a.imshow(niche_comp, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax_a.set_xticks(range(len(ct_names))); ax_a.set_xticklabels(ct_names, fontsize=4.5, rotation=90)
ax_a.set_yticks(range(K)); ax_a.set_yticklabels([f"N{n}" for n in range(K)], fontsize=6)
for ni in range(K):
    pct = ilc_niche_pct.get(ni, 0)
    top_ct = ct_names[np.argmax(niche_comp[ni])]
    ax_a.text(len(ct_names)+0.5, ni, f"{pct:.1f}% ILC", fontsize=5, color="#c44e52", va="center")
    ax_a.text(-1, ni, top_ct[:12], fontsize=4.5, color="#333", va="center", ha="right")
ax_a.set_title("a  Niche composition (log2 enrichment, % ILC-rich TLS)", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_a = plt.colorbar(im_a, ax=ax_a, shrink=0.7, pad=0.02)
cbar_a.set_label("log2(niche/global)", fontsize=6)

# ── Panel B: ILC-rich TLS distribution ──
ax_b = fig.add_axes([0.56, 0.68, 0.40, 0.28])
niches_sorted = sorted(ilc_niche_pct.keys(), key=lambda n: ilc_niche_pct[n], reverse=True)
pcts = [ilc_niche_pct[n] for n in niches_sorted]
colors_b = ["#c44e52" if p > 15 else "#4c72b0" if p > 8 else "#aaa" for p in pcts]
ax_b.bar(range(K), pcts, color=colors_b, alpha=0.7, edgecolor="black", linewidth=0.3)
ax_b.set_xticks(range(K)); ax_b.set_xticklabels([f"N{n}" for n in niches_sorted], fontsize=6)
ax_b.set_ylabel("% ILC-rich TLS spots", fontsize=7)
ax_b.axhline(y=100/K, color="#666", linestyle="--", linewidth=0.5, label=f"uniform ({100/K:.0f}%)")
ax_b.legend(fontsize=6, loc="upper right")
ax_b.set_title("b  ILC-rich TLS distribution across niches", fontsize=7.5, fontweight="bold", loc="left", pad=6)

# ── Panel C: Receptor detection ──
ax_c = fig.add_axes([0.06, 0.35, 0.42, 0.28])
df_rcp = df_all[df_all["type"]=="receptor"].nlargest(30, "pct_ILC").sort_values("pct_ILC")
rcp_mat = df_rcp[["pct_other","pct_ILC"]].values
im_c = ax_c.imshow(rcp_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_c.set_xticks([0,1]); ax_c.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_c.set_yticks(range(len(df_rcp)))
for i, (_, r) in enumerate(df_rcp.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax_c.text(-0.3, i, r["gene"], fontsize=5, fontstyle="italic", color=c, va="center", ha="right")
ax_c.set_title("c  Receptor detection rate (pct_expr)", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_c = plt.colorbar(im_c, ax=ax_c, shrink=0.7, pad=0.02)
cbar_c.set_label("pct_expr", fontsize=6)

# ── Panel D: Ligand detection ──
ax_d = fig.add_axes([0.56, 0.35, 0.40, 0.28])
df_lig = df_all[df_all["type"]=="ligand"].sort_values("pct_ILC")
lig_mat = df_lig[["pct_other","pct_ILC"]].values
im_d = ax_d.imshow(lig_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_d.set_xticks([0,1]); ax_d.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_d.set_yticks(range(len(df_lig)))
for i, (_, r) in enumerate(df_lig.iterrows()):
    ax_d.text(-0.3, i, r["gene"], fontsize=5, fontstyle="italic", color="#555", va="center", ha="right")
ax_d.set_title("d  Ligand detection rate", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_d = plt.colorbar(im_d, ax=ax_d, shrink=0.7, pad=0.02)
cbar_d.set_label("pct_expr", fontsize=6)

# ── Panel E: Receptor ranking ──
ax_e = fig.add_axes([0.06, 0.03, 0.90, 0.28])
df_rank = df_all[df_all["type"]=="receptor"].nlargest(30, "pct_ILC").sort_values("pct_ILC")
y_e = np.arange(len(df_rank))
colors_e = [CAT.get(r["category"], "#888") for _, r in df_rank.iterrows()]
ax_e.barh(y_e, df_rank["pct_ILC"].values, color=colors_e, alpha=0.6, height=0.7, edgecolor="black", linewidth=0.3)
for i, (_, r) in enumerate(df_rank.iterrows()):
    ax_e.text(r["pct_ILC"] + 0.005, i, r["gene"], fontsize=5.5, fontstyle="italic", va="center")
ax_e.set_yticks([])
ax_e.set_xlabel("ILC-rich TLS pct_expr", fontsize=7)
ax_e.set_xlim(0, 0.6)
ax_e.set_title("e  Receptor detection in ILC-rich TLS (top 30)", fontsize=7.5, fontweight="bold", loc="left", pad=6)
# Category legend
from matplotlib.patches import Patch
cats_in = set(df_rank["category"].values)
ax_e.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in],
            fontsize=5.5, loc="lower right", ncol=3, borderpad=0.3, title="Receptor category", title_fontsize=6)

# Save
fig.savefig(OUT / "fig_nature_5panel.svg", bbox_inches="tight")
fig.savefig(OUT / "fig_nature_5panel.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_nature_5panel.tiff", dpi=600, bbox_inches="tight")
plt.close()
print(f"Figure: {OUT}/fig_nature_5panel.{{svg,pdf,tiff}}")
for ni in range(K):
    print(f"  N{ni}: top={ct_names[np.argmax(niche_comp[ni])]}, ILC-rich={ilc_niche_pct[ni]:.1f}%")
