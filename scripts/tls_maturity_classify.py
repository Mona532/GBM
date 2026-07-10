"""TLS maturity classification following Pan-Cancer Atlas approach."""
import pandas as pd, numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"svg.fonttype":"none","pdf.fonttype":42,
    "font.size":8,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.6,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")
markers = pd.read_csv(r"E:/GBM/Pan-Cancer_Spatial_Atlas_TLS-main/data/TLS_signature/markers.tsv", sep="\t")
all_sig_genes = markers[markers["tls"]==1]["gene"].tolist()

meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")
nmf = pd.read_csv(ROOT / "tls_compnmf_rank5_unit_weights.csv")
spot_map = pd.read_csv(ROOT / "tls_component_spot_map.csv")

eco_names = {"E1":"Glial-CD4","E2":"TLS-structural","E3":"Vascular","E4":"Lymphocyte","E5":"Myeloid"}
eco_order = ["Lymphocyte","TLS-structural","Glial-CD4","Vascular","Myeloid"]
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00","#984ea3"]

import anndata as ad
scores = {}
for _, row in meta.iterrows():
    uid = row["unit_id"]; sid = row["sample"]
    h5_path = f"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated_tls16/{sid}.h5ad"
    try:
        a = ad.read_h5ad(h5_path)
        ge = a[:, a.var["feature_types"]=="Gene Expression"]
        spots = spot_map[spot_map["unit_id"]==uid]["barcode"].values
        spots = [s for s in spots if s in a.obs_names]
        if len(spots) < 5: continue
        ge_sub = ge[spots]
        X = ge_sub.X.toarray() if hasattr(ge_sub.X, "toarray") else ge_sub.X
        X_log = np.log1p(X)
        g2i = {g: i for i, g in enumerate(ge_sub.var_names)}
        gene_scores = [np.mean(X_log[:, g2i[g]]) for g in all_sig_genes if g in g2i]
        if gene_scores: scores[uid] = np.mean(gene_scores)
    except: pass

df = pd.DataFrame({"unit_id": list(scores.keys()), "tls_sig_score": list(scores.values())})
df["log_score"] = np.log1p(df["tls_sig_score"])
df = df.merge(nmf[["unit_id","dominant_ecotype"]], on="unit_id", how="left")
df["eco_name"] = df["dominant_ecotype"].map(eco_names)

# Cluster into 3 maturity levels
X = df[["log_score"]].values
X_scaled = StandardScaler().fit_transform(X)
labels = KMeans(n_clusters=3, random_state=42, n_init=10).fit_predict(X_scaled)
order = np.argsort([df[labels==i]["tls_sig_score"].mean() for i in range(3)])
label_map = {order[0]: "Immature", order[1]: "Primary", order[2]: "Secondary Mature"}
df["maturity"] = pd.Series(labels).map(label_map).values
df["maturity"] = pd.Categorical(df["maturity"], categories=["Immature","Primary","Secondary Mature"], ordered=True)

print("Maturity distribution:")
print(df["maturity"].value_counts().to_string())
print("\nPer-ecotype maturity (%):")
ct = pd.crosstab(df["eco_name"], df["maturity"])
ct = ct.loc[[e for e in eco_order if e in ct.index]]
print(ct.div(ct.sum(axis=1), axis=0).mul(100).round(1).to_string())

# Figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

# A: Boxplot
data = [df[df["eco_name"]==e]["tls_sig_score"].values for e in eco_order]
bp = ax1.boxplot(data, patch_artist=True, widths=0.5,
                  medianprops={"color":"black","linewidth":0.8}, flierprops={"marker":".","markersize":3})
for patch, c in zip(bp["boxes"], colors): patch.set_facecolor(c); patch.set_alpha(0.35)
ax1.set_xticklabels(eco_order, fontsize=9, rotation=25, ha="right")
ax1.set_ylabel("TLS signature score", fontsize=9)
ax1.set_title("TLS maturity score by ecotype", fontweight="bold", fontsize=10)

# B: Stacked bar (legend separated below)
ct_pct = ct.div(ct.sum(axis=1), axis=0)
ct_pct.plot(kind="bar", stacked=True, ax=ax2, color=["#2b83ba","#fdae61","#d7191c"], width=0.6, legend=False)
ax2.set_ylabel("Fraction", fontsize=9); ax2.set_xlabel("")
ax2.set_title("Maturity composition by ecotype", fontweight="bold", fontsize=10)
ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=9, rotation=25, ha="right")
ax2.set_ylim(0, 1.05)
# Legend below
handles = [plt.Rectangle((0,0),1,1,fc=c) for c in ["#2b83ba","#fdae61","#d7191c"]]
fig.legend(handles, ["Immature","Primary","Secondary Mature"], title="Maturity",
           fontsize=8, title_fontsize=9, loc="lower center", ncol=3, bbox_to_anchor=(0.65, 0.02))

fig.tight_layout(pad=2.0, rect=(0, 0.08, 1, 1))
fig.savefig(ROOT / "fig_tls_maturity_classify.jpg", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nSaved: {ROOT}/fig_tls_maturity_classify.jpg")
df.to_csv(ROOT / "tls_maturity_classification.csv", index=False)
