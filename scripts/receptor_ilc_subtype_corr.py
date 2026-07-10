"""ILC subtype stratification + Spearman correlation for receptor genes"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import spearmanr, wilcoxon, false_discovery_control
import warnings
warnings.filterwarnings("ignore")

SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
ILC_TYPES = ["ILC1","ILC2","ILC3"]
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna(): gene_cat[g] = cat_order[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]
OUT = Path(r"E:/GBM/results")

# Step 1: Global TLS P75 per ILC subtype
all_tls_ilc = {c: [] for c in ILC_TYPES}
for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() == 0: continue
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        for c in ILC_TYPES:
            all_tls_ilc[c].extend(q05[tls_mask, ct.index(c)].tolist())
GLOBAL_P75 = {c: np.percentile(all_tls_ilc[c], 75) for c in ILC_TYPES}
THRESH = {c: max(GLOBAL_P75[c], 1.0) for c in ILC_TYPES}
print(f"Thresholds: { {c: round(v,3) for c,v in THRESH.items()} }")

# Step 2: Per-ILC-subtype enrichment + Spearman correlation
subtype_deltas = []  # per-sample per-gene per-subtype log2FC
corr_rows = []       # Spearman within TLS spots

for tls_dir, h5_dir in DATASETS:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 6: continue
        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}

        # Gene expression: log-normalize
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr_raw = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values
        lib_size = expr_raw.sum(axis=1)
        expr_norm = np.log1p(expr_raw / (lib_size[:,None]/10000 + 1))

        non_expr = expr_norm[~tls_mask]

        # Define ILC-high per subtype
        for subtype in ILC_TYPES:
            si = ilc_idx[subtype]
            ilc_high = tls_mask & (q05[:, si] >= THRESH[subtype])
            if ilc_high.sum() < 3: continue
            ilc_expr = expr_norm[ilc_high]
            for g in ALL_RECEPTORS:
                if g not in vn: continue
                gidx = list(vn).index(g)
                tls_m = ilc_expr[:, gidx].mean()
                non_m = non_expr[:, gidx].mean()
                if tls_m == 0 or non_m == 0: continue
                l2fc = np.log2(tls_m / non_m)
                subtype_deltas.append({"sample": sd.name, "gene": g, "ILC_subtype": subtype,
                                        "log2FC": l2fc, "category": gene_cat.get(g,"?")})

        # Spearman correlation within TLS spots only
        tls_expr = expr_norm[tls_mask]
        tls_q05 = q05[tls_mask]
        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)
            ge_vals = tls_expr[:, gidx]
            if np.std(ge_vals) == 0: continue
            for c in ILC_TYPES:
                ilc_vals = tls_q05[:, ilc_idx[c]]
                if np.std(ilc_vals) == 0: continue
                rho, p = spearmanr(ge_vals, ilc_vals)
                corr_rows.append({"sample": sd.name, "gene": g, "ILC_subtype": c,
                                  "rho": rho, "p": p, "category": gene_cat.get(g,"?")})

# === Aggregate: per-subtype enrichment ===
df_sub = pd.DataFrame(subtype_deltas)
sub_agg = []
for g in sorted(df_sub["gene"].unique()):
    for st in ILC_TYPES:
        sub = df_sub[(df_sub["gene"]==g) & (df_sub["ILC_subtype"]==st)]
        if len(sub) < 5: continue
        l2fc = sub["log2FC"].dropna()
        if len(l2fc) < 5: continue
        _, pval = wilcoxon(l2fc, alternative="two-sided")
        sub_agg.append({"gene": g, "ILC_subtype": st, "n_samples": len(l2fc),
                         "median_log2FC": l2fc.median(), "pvalue": pval,
                         "category": gene_cat.get(g,"?")})
df_sub_agg = pd.DataFrame(sub_agg)
df_sub_agg["fdr"] = false_discovery_control(df_sub_agg["pvalue"].values)
df_sub_agg.to_csv(OUT / "receptor_ilc_subtype_enrichment.csv", index=False)

print(f"\n=== Per-subtype enrichment (top 5 per subtype, log2FC>1, FDR<0.1) ===")
for st in ILC_TYPES:
    sub = df_sub_agg[(df_sub_agg["ILC_subtype"]==st) & (df_sub_agg["median_log2FC"]>1) & (df_sub_agg["fdr"]<0.1)]
    sub = sub.nlargest(5, "median_log2FC")
    print(f"\n  {st}-high TLS (n_samples={sub['n_samples'].max()}):")
    for _, r in sub.iterrows():
        s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
        print(f"    {r['gene']:8s} {r['category']:22s} log2FC={r['median_log2FC']:+6.3f} FDR={r['fdr']:.4f} {s}")

# === Aggregate: Spearman correlation (per-gene, median rho across samples) ===
df_corr = pd.DataFrame(corr_rows)
corr_agg = []
for g in sorted(df_corr["gene"].unique()):
    for st in ILC_TYPES:
        sub = df_corr[(df_corr["gene"]==g) & (df_corr["ILC_subtype"]==st)]
        if len(sub) < 5: continue
        rhos = sub["rho"].dropna()
        if len(rhos) < 5: continue
        _, pval = wilcoxon(rhos, alternative="two-sided")
        corr_agg.append({"gene": g, "ILC_subtype": st, "n_samples": len(rhos),
                          "median_rho": rhos.median(), "pvalue": pval,
                          "category": gene_cat.get(g,"?")})
df_corr_agg = pd.DataFrame(corr_agg)
df_corr_agg["fdr"] = false_discovery_control(df_corr_agg["pvalue"].values)
df_corr_agg.to_csv(OUT / "receptor_spearman_ilc_subtype.csv", index=False)

print(f"\n=== Spearman correlation (top 10, |rho|>0.05, FDR<0.1) ===")
sig_corr = df_corr_agg[(df_corr_agg["fdr"]<0.1)].nlargest(15, "median_rho")
for _, r in sig_corr.iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} vs {r['ILC_subtype']:5s}  rho={r['median_rho']:+.4f}  FDR={r['fdr']:.4f}  n={r['n_samples']}  {s}")

print(f"\nOutput: {OUT}/receptor_ilc_subtype_enrichment.csv")
print(f"Output: {OUT}/receptor_spearman_ilc_subtype.csv")
