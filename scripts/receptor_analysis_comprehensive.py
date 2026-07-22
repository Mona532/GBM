"""Comprehensive receptor analysis: per-sample delta, log-norm, ILC subtype stratified, Spearman corr"""
import pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon, false_discovery_control
import warnings
warnings.filterwarnings("ignore")

SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
rx_df = pd.read_excel(r"E:\GBM\ti tianran2.xlsx")
cat_order = ["Excit (Glutamate)","Inhib (GABA/Gly)","Cholinergic (ACh)","DA/NE","Serotonin (5-HT)"]
gene_cat = {}
for idx, col in enumerate(rx_df.columns):
    for g in rx_df[col].dropna():
        gene_cat[g] = cat_order[idx]
ALL_RECEPTORS = sorted(gene_cat.keys())
ILC_TYPES = ["ILC1","ILC2","ILC3"]

DATASETS = [
    (Path(r"E:/GBM/results/tls_official_cut01"), Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_with_ilc")),
    (Path(r"E:/GBM/results/tls_visium_all"), Path(r"E:/GBM/ST_DATA/visium_all_h5ad")),
]

OUT = Path(r"E:/GBM/results")

# ====== Step 1: Long table of per-spot expression ======
long_rows = []
sample_deltas = []  # per-sample per-gene delta_pct, delta_mean

for ds_name, tls_dir, h5_dir in [("D1", *DATASETS[0]), ("D2", *DATASETS[1])]:
    for sd in sorted(tls_dir.iterdir()):
        if not sd.is_dir(): continue
        if sd.name in SKIP_DMG: continue
        tls_csv = sd / "tls_spot_scores_official_relaxed.csv"
        h5 = h5_dir / f"{sd.name}.h5ad"
        if not tls_csv.exists() or not h5.exists(): continue
        tls = pd.read_csv(tls_csv)
        tls_mask = (tls["TLS.region"]=="TLS").values
        if tls_mask.sum() < 3: continue

        adata = ad.read_h5ad(h5)
        q05 = adata.obsm["c2l_ilc_q05"]
        if hasattr(q05,"values"): q05 = q05.values
        ct = list(adata.uns["c2l_ilc_cell_types"])

        # ILC subtypes
        ilc_idx = {c: ct.index(c) for c in ILC_TYPES}
        ilc_vals = {c: q05[:, ilc_idx[c]] for c in ILC_TYPES}

        # Gene expression: log-normalize
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"] if "feature_types" in adata.var else adata
        expr_raw = ge.X.toarray() if hasattr(ge.X,"toarray") else ge.X
        vn = ge.var_names.values
        # Simple log-normalization: log1p(CP10k)
        lib_size = expr_raw.sum(axis=1)
        expr_norm = np.log1p(expr_raw / (lib_size[:,None] / 10000 + 1))

        # Groups
        tls_group = np.where(tls_mask, "TLS", "nonTLS")

        for g in ALL_RECEPTORS:
            if g not in vn: continue
            gidx = list(vn).index(g)

            # Per-spot rows
            for i in range(len(tls)):
                long_rows.append({
                    "spot": adata.obs_names[i], "sample": sd.name, "dataset": ds_name,
                    "group": tls_group[i],
                    "ILC1_q05": ilc_vals["ILC1"][i], "ILC2_q05": ilc_vals["ILC2"][i], "ILC3_q05": ilc_vals["ILC3"][i],
                    "gene": g, "category": gene_cat.get(g,"?"),
                    "expr": expr_norm[i, gidx],
                })

            # Per-sample delta (TLS vs nonTLS)
            tls_expr = expr_norm[tls_mask, gidx]
            non_expr = expr_norm[~tls_mask, gidx]
            if len(tls_expr) < 3 or len(non_expr) < 3: continue
            sample_deltas.append({
                "sample": sd.name, "dataset": ds_name, "gene": g,
                "pct_TLS": (tls_expr > 0).mean(), "pct_nonTLS": (non_expr > 0).mean(),
                "mean_TLS": tls_expr.mean(), "mean_nonTLS": non_expr.mean(),
                "delta_pct": (tls_expr > 0).mean() - (non_expr > 0).mean(),
                "delta_mean": tls_expr.mean() - non_expr.mean(),
            })

df_long = pd.DataFrame(long_rows)
df_delta = pd.DataFrame(sample_deltas)
df_long.to_csv(OUT / "receptor_long_table.csv", index=False)
df_delta.to_csv(OUT / "receptor_sample_deltas.csv", index=False)
print(f"Long table: {len(df_long)} rows, {df_long['sample'].nunique()} samples, {df_long['gene'].nunique()} genes")
print(f"Per-sample deltas: {len(df_delta)} rows, {df_delta['sample'].nunique()} samples")

# ====== Step 2: Aggregate per gene (per-sample delta → Wilcoxon test) ======
agg_rows = []
for g in sorted(df_delta["gene"].unique()):
    sub = df_delta[df_delta["gene"] == g]
    if len(sub) < 5: continue
    # Test delta_pct > 0
    _, p_pct = wilcoxon(sub["delta_pct"].values, alternative="two-sided")
    _, p_mean = wilcoxon(sub["delta_mean"].values, alternative="two-sided")
    agg_rows.append({
        "gene": g, "n_samples": len(sub),
        "median_delta_pct": sub["delta_pct"].median(),
        "mean_delta_pct": sub["delta_pct"].mean(),
        "pct_TLS_avg": sub["pct_TLS"].mean(),
        "pct_nonTLS_avg": sub["pct_nonTLS"].mean(),
        "median_delta_mean": sub["delta_mean"].median(),
        "p_pct": p_pct, "p_mean": p_mean,
        "category": gene_cat.get(g, "?"),
    })

df_agg = pd.DataFrame(agg_rows)
df_agg["fdr_pct"] = false_discovery_control(df_agg["p_pct"].values)
df_agg["fdr_mean"] = false_discovery_control(df_agg["p_mean"].values)
df_agg = df_agg.sort_values("median_delta_pct", ascending=False)
df_agg.to_csv(OUT / "receptor_aggregated_stats.csv", index=False)

print(f"\n=== Aggregate per-gene stats ===")
sig_pct = df_agg[df_agg["fdr_pct"] < 0.1]
sig_mean = df_agg[df_agg["fdr_mean"] < 0.1]
print(f"Significant delta_pct (FDR<0.1): {len(sig_pct)}")
print(f"Significant delta_mean (FDR<0.1): {len(sig_mean)}")

print("\nTop by delta_pct (FDR<0.1):")
for _, r in sig_pct.head(10).iterrows():
    s = "***" if r["fdr_pct"]<0.001 else "**" if r["fdr_pct"]<0.01 else "*"
    print(f"  {r['gene']:8s} {r['category']:22s}  Δpct={r['median_delta_pct']:+.4f}  TLS={r['pct_TLS_avg']:.3f}  nonTLS={r['pct_nonTLS_avg']:.3f}  FDR={r['fdr_pct']:.4f}  n={r['n_samples']}  {s}")

# ====== Step 3: Spearman correlation within TLS spots ======
corr_rows = []
for g in sorted(df_long["gene"].unique()):
    sub_tls = df_long[df_long["group"] == "TLS"]
    gene_tls = sub_tls[sub_tls["gene"] == g]
    if len(gene_tls) < 50: continue
    for c in ILC_TYPES:
        rho, p = spearmanr(gene_tls["expr"], gene_tls[f"{c}_q05"])
        corr_rows.append({"gene": g, "ILC_type": c, "rho": rho, "p": p, "n": len(gene_tls)})

df_corr = pd.DataFrame(corr_rows)
df_corr["fdr"] = false_discovery_control(df_corr["p"].values)

# Per gene: strongest ILC
gene_strongest = df_corr.groupby("gene").apply(lambda x: x.loc[x["rho"].idxmax(), "ILC_type"]).reset_index(name="strongest_ILC")
df_corr = df_corr.merge(gene_strongest, on="gene")
df_corr.to_csv(OUT / "receptor_spearman_corr.csv", index=False)

print(f"\n=== Spearman correlation (top by |rho|, FDR<0.1) ===")
sig_corr = df_corr[df_corr["fdr"] < 0.1].sort_values("rho", ascending=False)
for _, r in sig_corr.head(15).iterrows():
    s = "***" if r["fdr"]<0.001 else "**" if r["fdr"]<0.01 else "*"
    print(f"  {r['gene']:8s} vs {r['ILC_type']:5s}  rho={r['rho']:+.4f}  FDR={r['fdr']:.4f}  n={r['n']}  strongest={r['strongest_ILC']}  {s}")

print(f"\nAll output in: {OUT}")
