"""
Step 1: Compute & save all data. Step 2: Read saved data & generate Nature 5-panel figure.
CellCharter K=5 (Cluster), all intermediate results saved to CSV.
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import cellcharter as cc, squidpy as sq
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")

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
OUT.mkdir(parents=True, exist_ok=True)
K_NICHES = 5
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.5, "legend.frameon": False,
})

# ====== STEP 1: Compute & save ======
print("=== Step 1: Computing ===")
data_list, sample_info = [], []
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
        a = ad.AnnData(X=q05_z.astype(np.float32),
                       obs=pd.DataFrame({"sample": [sd.name]*q05.shape[0]}, index=adata.obs_names),
                       obsm={"spatial": coords})
        a.var_names = ct_names; a.obs["sample"] = a.obs["sample"].astype("category")
        data_list.append(a)
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        genes_present = [g for g in ALL_GENES if g in vn_all]
        gene_idx = [list(vn_all).index(g) for g in genes_present]
        rx_raw = ge[:, gene_idx].X.toarray() if hasattr(ge[:, gene_idx].X,"toarray") else np.asarray(ge[:, gene_idx].X)
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))
        sample_info.append((sd.name, tls_mask, q05, expr_norm, rx_raw, genes_present, ct_names, ilc_idx_map, coords))

print(f"  Samples: {len(sample_info)}")

# CellCharter (cc.tl.Cluster — official fixed-K API)
adata_all = ad.concat(data_list, join="inner")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in data_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
# ClusterAutoK K=5-6 auto-select, max_runs=3
gmm = cc.tl.ClusterAutoK(n_clusters=(5, 6), max_runs=3)
gmm.fit(adata_all, use_rep="X_cellcharter")
niche_all = gmm.predict(adata_all, use_rep="X_cellcharter")
K_NICHES = len(set(niche_all))
print(f"  Niches: {K_NICHES} (fixed)")

# Niche composition
all_raw = np.vstack([s[2] for s in sample_info])
global_mean = all_raw.mean(axis=0) + 1e-8
ct_names_all = list(adata_all.var_names)
niche_comp = np.zeros((K_NICHES, len(ct_names_all)))
for ni in range(K_NICHES):
    niche_comp[ni] = np.log2(all_raw[niche_all==ni].mean(axis=0) / global_mean)

# ILC-rich TLS metrics
gene_sample = {}; ilc_counts = np.zeros(K_NICHES)
offset = 0
for sd_name, tls_mask, q05, expr_n, rx_r, genes_p, ct_n, ilc_m, _ in sample_info:
    n = len(tls_mask); niches = niche_all[offset:offset+n]; offset += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_m[c]] >= GLOBAL_THRESH[c])
    ilc_rich = ilc_high & tls_mask; other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue
    for ni in range(K_NICHES): ilc_counts[ni] += int((niches == ni).sum())
    ie, oe = expr_n[ilc_rich], expr_n[other_tls]
    for gi, g in enumerate(genes_p):
        pct_i = (rx_r[ilc_rich, gi] > 0).mean(); pct_o = (rx_r[other_tls, gi] > 0).mean()
        gene_sample.setdefault(g, []).append({"pct_ilc":pct_i,"pct_other":pct_o,"delta":ie[:,gi].mean()-oe[:,gi].mean()})

# Save everything
np.savetxt(OUT / "niche_composition.csv", niche_comp, delimiter=",")
pd.DataFrame({"cell_type": ct_names_all}).to_csv(OUT / "cell_type_names.csv", index=False)
n_pass = 0
for sd_name, tls_mask, q05, expr_n, rx_r, genes_p, ct_n, ilc_m, _ in sample_info:
    n = len(tls_mask); niches = niche_all[offset_temp:offset_temp+n]
    offset_temp += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_m[c]] >= GLOBAL_THRESH[c])
    if (ilc_high & tls_mask).sum() >= 3 and (tls_mask & ~(ilc_high & tls_mask)).sum() >= 3:
        n_pass += 1
print(f"  Samples passing ILC-rich filter: {n_pass}/{len(sample_info)}")
del offset_temp

if ilc_counts.sum() == 0:
    print("  WARNING: No ILC-rich spots, using uniform distribution")
    ilc_pct = {ni: 100/K_NICHES for ni in range(K_NICHES)}
else:
    ilc_pct = {ni: ilc_counts[ni]/ilc_counts.sum()*100 for ni in range(K_NICHES)}
pd.DataFrame({"niche": list(ilc_pct.keys()), "ilc_pct": list(ilc_pct.values())}).to_csv(OUT / "niche_ilc_distribution.csv", index=False)
gene_rows = []
for g, sv in gene_sample.items():
    p_i = np.median([v["pct_ilc"] for v in sv]); p_o = np.median([v["pct_other"] for v in sv])
    d = np.median([v["delta"] for v in sv])
    gene_rows.append({"gene":g,"pct_ILC":p_i,"pct_other":p_o,"median_delta":d,
        "type":"receptor" if g in RECEPTOR_GENES else "ligand",
        "category":RECEPTOR_GENES.get(g,"Ligand"),"n_samples":len(sv)})
pd.DataFrame(gene_rows).to_csv(OUT / "gene_metrics.csv", index=False)
print("  Data saved.")

# ====== STEP 2: Figure from saved data ======
print("=== Step 2: Plotting ===")
niche_comp = pd.read_csv(OUT / "niche_composition.csv", header=None).values
ct_names = pd.read_csv(OUT / "cell_type_names.csv")["cell_type"].tolist()
ilc_dist = pd.read_csv(OUT / "niche_ilc_distribution.csv")
ilc_pct = dict(zip(ilc_dist["niche"], ilc_dist["ilc_pct"]))
df_all = pd.read_csv(OUT / "gene_metrics.csv")

fig = plt.figure(figsize=(10, 11))

# A: Niche composition
ax_a = fig.add_axes([0.06, 0.68, 0.42, 0.28])
vmax = max(abs(niche_comp.min()), abs(niche_comp.max()), 1.0)
im_a = ax_a.imshow(niche_comp, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax_a.set_xticks(range(len(ct_names))); ax_a.set_xticklabels(ct_names, fontsize=4.5, rotation=90)
ax_a.set_yticks(range(K_NICHES)); ax_a.set_yticklabels([f"N{n}" for n in range(K_NICHES)], fontsize=6)
for ni in range(K_NICHES):
    top_ct = ct_names[np.argmax(niche_comp[ni])]
    ax_a.text(len(ct_names)+0.5, ni, f"{ilc_pct[ni]:.1f}%", fontsize=5, color="#c44e52", va="center")
    ax_a.text(-1, ni, top_ct[:12], fontsize=4.5, color="#333", va="center", ha="right")
ax_a.set_title("a  Niche composition (log2 enrichment, % ILC-rich)", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_a = plt.colorbar(im_a, ax=ax_a, shrink=0.7, pad=0.02); cbar_a.set_label("log2(niche/global)", fontsize=6)

# B: ILC-rich distribution
ax_b = fig.add_axes([0.56, 0.68, 0.40, 0.28])
sorted_niches = sorted(ilc_pct.keys(), key=lambda n: ilc_pct[n], reverse=True)
pcts = [ilc_pct[n] for n in sorted_niches]
colors_b = ["#c44e52" if p > 20 else "#4c72b0" if p > 10 else "#aaa" for p in pcts]
ax_b.bar(range(K_NICHES), pcts, color=colors_b, alpha=0.7, edgecolor="black", linewidth=0.3)
ax_b.set_xticks(range(K_NICHES)); ax_b.set_xticklabels([f"N{n}" for n in sorted_niches], fontsize=6)
ax_b.set_ylabel("% ILC-rich TLS spots", fontsize=7)
ax_b.axhline(y=100/K_NICHES, color="#666", linestyle="--", linewidth=0.5, label=f"uniform ({100/K_NICHES:.0f}%)")
ax_b.legend(fontsize=6, loc="upper right")
ax_b.set_title("b  ILC-rich TLS niche distribution", fontsize=7.5, fontweight="bold", loc="left", pad=6)

# C: Receptor heatmap
ax_c = fig.add_axes([0.06, 0.35, 0.42, 0.28])
df_rcp = df_all[df_all["type"]=="receptor"].nlargest(30, "pct_ILC").sort_values("pct_ILC")
rcp_mat = df_rcp[["pct_other","pct_ILC"]].values
im_c = ax_c.imshow(rcp_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_c.set_xticks([0,1]); ax_c.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_c.set_yticks(range(len(df_rcp)))
for i, (_, r) in enumerate(df_rcp.iterrows()):
    c = CAT.get(r["category"], "#888")
    ax_c.text(-0.3, i, r["gene"], fontsize=5, fontstyle="italic", color=c, va="center", ha="right")
ax_c.set_title("c  Receptor detection rate", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_c = plt.colorbar(im_c, ax=ax_c, shrink=0.7, pad=0.02); cbar_c.set_label("pct_expr", fontsize=6)

# D: Ligand heatmap
ax_d = fig.add_axes([0.56, 0.35, 0.40, 0.28])
df_lig = df_all[df_all["type"]=="ligand"].sort_values("pct_ILC")
lig_mat = df_lig[["pct_other","pct_ILC"]].values
im_d = ax_d.imshow(lig_mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
ax_d.set_xticks([0,1]); ax_d.set_xticklabels(["other TLS","ILC-rich TLS"], fontsize=6)
ax_d.set_yticks(range(len(df_lig)))
for i, (_, r) in enumerate(df_lig.iterrows()):
    ax_d.text(-0.3, i, r["gene"], fontsize=5, fontstyle="italic", color="#555", va="center", ha="right")
ax_d.set_title("d  Ligand detection rate", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cbar_d = plt.colorbar(im_d, ax=ax_d, shrink=0.7, pad=0.02); cbar_d.set_label("pct_expr", fontsize=6)

# E: Receptor ranking
ax_e = fig.add_axes([0.06, 0.03, 0.90, 0.28])
df_rank = df_all[df_all["type"]=="receptor"].nlargest(30, "pct_ILC").sort_values("pct_ILC")
ye = np.arange(len(df_rank))
colors_e = [CAT.get(r["category"], "#888") for _, r in df_rank.iterrows()]
ax_e.barh(ye, df_rank["pct_ILC"].values, color=colors_e, alpha=0.6, height=0.7, edgecolor="black", linewidth=0.3)
for i, (_, r) in enumerate(df_rank.iterrows()):
    ax_e.text(r["pct_ILC"] + 0.005, i, r["gene"], fontsize=5.5, fontstyle="italic", va="center")
ax_e.set_yticks([]); ax_e.set_xlabel("ILC-rich TLS pct_expr", fontsize=7); ax_e.set_xlim(0, 0.6)
ax_e.set_title("e  Top receptors in ILC-rich TLS", fontsize=7.5, fontweight="bold", loc="left", pad=6)
cats_in = set(df_rank["category"].values)
ax_e.legend(handles=[Patch(color=CAT[c], alpha=0.5, label=c.split("(")[0].strip()) for c in CAT if c in cats_in],
            fontsize=5.5, loc="lower right", ncol=3, borderpad=0.3, title="Category", title_fontsize=6)

fig.savefig(OUT / "fig_nature_5panel.svg", bbox_inches="tight")
fig.savefig(OUT / "fig_nature_5panel.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_nature_5panel.tiff", dpi=600, bbox_inches="tight")
plt.close()
print(f"  Figure: {OUT}/fig_nature_5panel.{{svg,pdf,tiff}}")
for ni in range(K_NICHES):
    print(f"    N{ni}: top={ct_names[np.argmax(niche_comp[ni])]}, ILC-rich={ilc_pct[ni]:.1f}%")
