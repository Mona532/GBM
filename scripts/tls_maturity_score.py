"""TLS maturity scoring using Pan-Cancer Atlas signature genes."""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"svg.fonttype":"none","pdf.fonttype":42,
    "font.size":8,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.6,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

# Load paper's TLS signature
markers = pd.read_csv(r"E:/GBM/Pan-Cancer_Spatial_Atlas_TLS-main/data/TLS_signature/markers.tsv", sep="\t")
tls_genes = markers[markers["tls"]==1]["gene"].tolist()
lymphoid_genes = markers[markers["lymphoid"]==1]["gene"].tolist()

# Load pseudobulk counts
counts = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")  # need gene counts from R
# Actually load from C2L - use per-component c2l spot map to aggregate
spot_map = pd.read_csv(ROOT / "tls_component_spot_map.csv")
nmf = pd.read_csv(ROOT / "tls_compnmf_rank5_unit_weights.csv")
meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")

# Load gene expression from one sample's h5ad to check gene availability
import anndata as ad
sample_h5 = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated_tls16/AT10-BRA-5-FO-1_1.h5ad"
a = ad.read_h5ad(sample_h5)
ge_mask = a.var["feature_types"] == "Gene Expression"
spatial_genes = a.var_names[ge_mask].tolist()

tls_avail = [g for g in tls_genes if g in spatial_genes]
lym_avail = [g for g in lymphoid_genes if g in spatial_genes]
print(f"TLS genes in spatial: {len(tls_avail)}/{len(tls_genes)}")
print(f"Lymphoid genes in spatial: {len(lym_avail)}/{len(lymphoid_genes)}")

# Score each component: sum of normalized gene expression across all spots
# For each sample, compute per-spot log1p expression, then component-level mean
scores = []
for uid in meta["unit_id"]:
    sid = meta[meta["unit_id"]==uid]["sample"].values[0]
    h5_path = f"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated_tls16/{sid}.h5ad"
    try:
        a2 = ad.read_h5ad(h5_path)
        ge2 = a2[:, a2.var["feature_types"]=="Gene Expression"]
        # Get spots for this component
        spots = spot_map[spot_map["unit_id"]==uid]["barcode"].values
        spots = [s for s in spots if s in a2.obs_names]
        if len(spots) < 5: continue
        ge_sub = ge2[spots]
        # Log-normalize
        X = ge_sub.X.toarray() if hasattr(ge_sub.X, "toarray") else ge_sub.X
        X_log = np.log1p(X)
        gene_idx = {g: i for i, g in enumerate(ge_sub.var_names)}
        tls_score = np.mean([X_log[:, gene_idx[g]].mean() for g in tls_avail if g in gene_idx])
        lym_score = np.mean([X_log[:, gene_idx[g]].mean() for g in lym_avail if g in gene_idx])
        scores.append({"unit_id": uid, "tls_score": tls_score, "lymphoid_score": lym_score, "n_spots": len(spots)})
    except Exception as e:
        pass

df = pd.DataFrame(scores)
df.to_csv(ROOT / "tls_maturity_scores.csv", index=False)
print(f"Scored {len(df)} components")
print(df.describe())

# Merge with ecotype
df_eco = df.merge(nmf[["unit_id","dominant_ecotype"]], on="unit_id", how="left")
print("\nPer-ecotype maturity scores:")
print(df_eco.groupby("dominant_ecotype")[["tls_score","lymphoid_score"]].mean().round(3).to_string())

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
eco_order = sorted(df_eco["dominant_ecotype"].dropna().unique())
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00","#984ea3"]
for i, e in enumerate(eco_order):
    sub = df_eco[df_eco["dominant_ecotype"]==e]
    ax1.scatter(sub["tls_score"], sub["lymphoid_score"], s=8, alpha=0.4, color=colors[i], label=f"E{e} (n={len(sub)})")
ax1.set_xlabel("TLS signature score"); ax1.set_ylabel("Lymphoid signature score")
ax1.legend(fontsize=7); ax1.set_title("TLS maturity by ecotype", fontweight="bold")

# Boxplot of TLS score per ecotype
data = [df_eco[df_eco["dominant_ecotype"]==e]["tls_score"].values for e in eco_order]
bp = ax2.boxplot(data, patch_artist=True, widths=0.5)
for patch, c in zip(bp["boxes"], colors): patch.set_facecolor(c); patch.set_alpha(0.3)
ax2.set_xticklabels(eco_order, fontsize=8)
ax2.set_ylabel("TLS signature score")
ax2.set_title("TLS signature by ecotype", fontweight="bold")
fig.tight_layout()
fig.savefig(ROOT / "fig_tls_maturity.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nSaved: {ROOT}/fig_tls_maturity.jpg")
