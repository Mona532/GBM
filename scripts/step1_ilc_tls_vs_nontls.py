"""
ILC-high TLS vs ILC-high non-TLS (downsampled) — receptor/ligand comparison
P90 thresholds, per-sample balance by downsampling non-TLS group
"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
CAT_LIST = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
RECEPTOR_GENES = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): RECEPTOR_GENES[g] = CAT_LIST[idx]
LIGANDS = ["NMU","VIP","ADM","CALCA","CALCB","NPY","TAC1","TAC3","NTS","CCK","GRP","GAL","PENK","PDYN","PNOC","CRH","UCN","AGRP","POMC","MCH"]
ALL_GENES = sorted(set(list(RECEPTOR_GENES.keys()) + LIGANDS))
ILC_TYPES = ["ILC1","ILC2","ILC3"]
P90 = {"ILC1":0.707,"ILC2":0.709,"ILC3":0.753}
SKIP = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results"); OUT.mkdir(parents=True, exist_ok=True)
DS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

# Per-sample: downsample nonTLS ILC-high to match TLS ILC-high count, compute per-gene metrics
gene_sample = {}
for tls_dir, h5_dir in DS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir() or sd.name in SKIP: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"; h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        if "barcode" in tls.columns: tls = tls.set_index("barcode")
        adata = ad.read_h5ad(h5)
        shared = adata.obs_names.intersection(tls.index)
        if len(shared) < 100: continue
        adata = adata[shared]; tls = tls.loc[shared]
        tls_mask = (tls["TLS.region"]=="TLS").values
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}
        ilc_high = np.zeros(len(tls), dtype=bool)
        for c in ILC_TYPES: ilc_high |= (q05[:, ilc_idx[c]] >= P90[c])
        tls_ilc = ilc_high & tls_mask
        nontls_ilc = ilc_high & ~tls_mask
        n_tls = tls_ilc.sum(); n_non = nontls_ilc.sum()
        if n_tls < 3 or n_non < 3: continue

        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        lib_size = np.array(ge.X.sum(axis=1)).flatten() + 1
        vn_all = ge.var_names.values
        genes_present = [g for g in ALL_GENES if g in vn_all]
        gene_idx = [list(vn_all).index(g) for g in genes_present]
        rx_raw = ge[:, gene_idx].X.toarray() if hasattr(ge[:, gene_idx].X,"toarray") else np.asarray(ge[:, gene_idx].X)
        expr_norm = np.log1p(rx_raw / (lib_size[:,None]/10000))

        # Downsample nonTLS to match TLS count
        rng = np.random.RandomState(42)
        nontls_idx = np.where(nontls_ilc)[0]
        n_sample = min(n_tls, n_non)
        nontls_sampled = rng.choice(nontls_idx, n_sample, replace=False)
        tls_idx = np.where(tls_ilc)[0]

        ie = expr_norm[tls_idx]; ne = expr_norm[nontls_sampled]
        for gi, g in enumerate(genes_present):
            pct_t = (rx_raw[tls_idx, gi] > 0).mean(); pct_n = (rx_raw[nontls_sampled, gi] > 0).mean()
            m_t = ie[:, gi].mean(); m_n = ne[:, gi].mean()
            gene_sample.setdefault(g, []).append({"pct_TLS":pct_t,"pct_nonTLS":pct_n,"mean_TLS":m_t,"mean_nonTLS":m_n,"delta":m_t-m_n})

# Aggregate
rows = []
for g, sv in gene_sample.items():
    n = len(sv)
    pct_t = np.array([v["pct_TLS"] for v in sv]); pct_n = np.array([v["pct_nonTLS"] for v in sv])
    m_t = np.array([v["mean_TLS"] for v in sv]); m_n = np.array([v["mean_nonTLS"] for v in sv])
    d = m_t - m_n
    mask_det = pct_t > 0
    rows.append({"gene":g,"n_total":n,"n_detected":int(mask_det.sum()),
        "detect_rate":mask_det.mean(),
        "pct_TLS_detected":np.median(pct_t[mask_det]) if mask_det.sum()>0 else 0,
        "pct_nonTLS_detected":np.median(pct_n[mask_det]) if mask_det.sum()>0 else 0,
        "mean_TLS_detected":np.median(m_t[mask_det]) if mask_det.sum()>0 else 0,
        "mean_nonTLS_detected":np.median(m_n[mask_det]) if mask_det.sum()>0 else 0,
        "median_delta_all":np.median(d),"prop_enriched":(d>0).mean(),
        "type":"receptor" if g in RECEPTOR_GENES else "ligand",
        "category":RECEPTOR_GENES.get(g,"Ligand")})

df = pd.DataFrame(rows)
df.to_csv(OUT/"gene_metrics_ilc_tls_vs_nontls.csv", index=False)
print(f"Genes: {len(df)}, receptors: {(df['type']=='receptor').sum()}, ligands: {(df['type']=='ligand').sum()}")
dr = df[df['type']=='receptor']
print(f"Receptors with detect_rate>=0.1: {(dr['detect_rate']>=0.1).sum()}, prop+>=0.6: {(dr['prop_enriched']>=0.6).sum()}")
cand = dr[(dr['detect_rate']>=0.1)&(dr['prop_enriched']>=0.5)].sort_values('median_delta_all',ascending=False)
print(f"Top candidates:")
for _,r in cand.head(8).iterrows():
    print(f"  {r['gene']:8s} detect={r['detect_rate']:.2f} prop+={r['prop_enriched']:.2f} delta={r['median_delta_all']:+.4f}")
