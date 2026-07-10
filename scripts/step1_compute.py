"""Step 1: CellCharter K=5/6 + ILC-rich TLS metrics → save all CSVs"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import cellcharter as cc, squidpy as sq
import warnings; warnings.filterwarnings("ignore")

rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
CAT_LIST = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
RECEPTOR_GENES = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): RECEPTOR_GENES[g] = CAT_LIST[idx]
LIGANDS = ["NMU","VIP","ADM","CALCA","CALCB","NPY","TAC1","TAC3","NTS","CCK","GRP","GAL","PENK","PDYN","PNOC","CRH","UCN","AGRP","POMC","MCH"]
ALL_GENES = sorted(set(list(RECEPTOR_GENES.keys()) + LIGANDS))
ILC_TYPES = ["ILC1","ILC2","ILC3"]
THRESH = {"ILC1":1.034,"ILC2":1.0,"ILC3":1.035}
SKIP = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results"); OUT.mkdir(parents=True, exist_ok=True)
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

print("Loading samples...")
data_list, sample_info = [], []
for tls_dir, h5_dir in DS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir() or sd.name in SKIP: continue
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
        ilc_idx = {c: ct_names.index(c) for c in ILC_TYPES}
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
        sample_info.append((sd.name, tls_mask, q05, expr_norm, rx_raw, genes_present, ct_names, ilc_idx, coords))
print(f"  {len(sample_info)} samples")

# CellCharter
adata_all = ad.concat(data_list, join="inner")
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")
adata_all.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in data_list])
sq.gr.spatial_neighbors(adata_all, library_key="sample", coord_type="generic", delaunay=True)
cc.gr.aggregate_neighbors(adata_all, n_layers=1, use_rep=None, out_key="X_cellcharter", sample_key="sample")
gmm = cc.tl.ClusterAutoK(n_clusters=(5, 6), max_runs=3)
gmm.fit(adata_all, use_rep="X_cellcharter")
niche_all = gmm.predict(adata_all, use_rep="X_cellcharter").astype(int)
K = len(set(niche_all))
print(f"  Niches: {K}, dtype={niche_all.dtype}, unique={sorted(set(niche_all))[:8]}")

# Niche composition
all_raw = np.vstack([s[2] for s in sample_info])
global_mean = all_raw.mean(axis=0) + 1e-8
ct_names_all = list(adata_all.var_names)
niche_comp = np.zeros((K, len(ct_names_all)))
for ni in range(K):
    niche_comp[ni] = np.log2(all_raw[niche_all==ni].mean(axis=0) / global_mean)

# ILC-rich TLS + gene metrics
gene_sample = {}; ilc_counts = np.zeros(K); n_pass = 0
offset = 0
for sd_name, tls_mask, q05, expr_n, rx_r, genes_p, ct_n, ilc_m, _ in sample_info:
    n = len(tls_mask)
    niches = niche_all[offset:offset+n]; offset += n
    ilc_high = np.zeros(n, dtype=bool)
    for c in ILC_TYPES: ilc_high |= (q05[:, ilc_m[c]] >= THRESH[c])
    ilc_rich = ilc_high & tls_mask
    other_tls = tls_mask & ~ilc_rich
    if ilc_rich.sum() < 3 or other_tls.sum() < 3: continue
    n_pass += 1
    for ni in range(K): ilc_counts[ni] += int((niches == ni).sum())
    ie, oe = expr_n[ilc_rich], expr_n[other_tls]
    for gi, g in enumerate(genes_p):
        pct_i = (rx_r[ilc_rich, gi] > 0).mean()
        pct_o = (rx_r[other_tls, gi] > 0).mean()
        m_i = ie[:, gi].mean()
        m_o = oe[:, gi].mean()
        gene_sample.setdefault(g, []).append({"pct_ilc":pct_i,"pct_other":pct_o,"delta":m_i-m_o,"mean_ILC":m_i,"mean_other":m_o})
print(f"  ILC-rich filter passed: {n_pass}/{len(sample_info)}, ilc_counts={ilc_counts}")

# Save
np.savetxt(OUT / "niche_composition.csv", niche_comp, delimiter=",")
pd.DataFrame({"cell_type": ct_names_all}).to_csv(OUT / "cell_type_names.csv", index=False)
if ilc_counts.sum() > 0:
    ilc_pct = {ni: ilc_counts[ni]/ilc_counts.sum()*100 for ni in range(K)}
else:
    ilc_pct = {ni: 100/K for ni in range(K)}
pd.DataFrame({"niche": list(ilc_pct.keys()), "ilc_pct": list(ilc_pct.values())}).to_csv(OUT / "niche_ilc_distribution.csv", index=False)

gene_rows = []
for g, sv in gene_sample.items():
    n_total = len(sv)
    pct_i_all = np.array([v["pct_ilc"] for v in sv])
    pct_o_all = np.array([v["pct_other"] for v in sv])
    m_i_all = np.array([v["mean_ILC"] for v in sv])
    m_o_all = np.array([v["mean_other"] for v in sv])
    d_all = m_i_all - m_o_all
    # Detection
    n_detected = int((pct_i_all > 0).sum())
    detect_rate = n_detected / n_total
    # Among detected samples
    mask_det = pct_i_all > 0
    pct_i_det = np.median(pct_i_all[mask_det]) if mask_det.sum() > 0 else 0
    m_i_det = np.median(m_i_all[mask_det]) if mask_det.sum() > 0 else 0
    d_det = np.median(d_all[mask_det]) if mask_det.sum() > 0 else 0
    # Among all samples
    d_all_med = np.median(d_all)
    # Enrichment direction
    prop_enriched = (d_all > 0).mean()
    pct_o_det = np.median(pct_o_all[mask_det]) if mask_det.sum() > 0 else 0
    m_o_det = np.median(m_o_all[mask_det]) if mask_det.sum() > 0 else 0
    gene_rows.append({"gene":g, "n_total":n_total, "n_detected":n_detected,
        "detect_rate":detect_rate, "pct_ILC_detected":pct_i_det,
        "pct_other_detected":pct_o_det, "mean_ILC_detected":m_i_det,
        "mean_other_detected":m_o_det, "median_delta_detected":d_det,
        "median_delta_all":d_all_med, "prop_enriched":prop_enriched,
        "type":"receptor" if g in RECEPTOR_GENES else "ligand",
        "category":RECEPTOR_GENES.get(g,"Ligand")})
pd.DataFrame(gene_rows).to_csv(OUT / "gene_metrics.csv", index=False)
print("  All CSVs saved.")
