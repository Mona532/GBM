"""
All final figures — using merged c2l + NMF data (post CD4_T + Plasma fix).
"""
import pandas as pd, numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":8,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.6,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

# ====== Load data ======
basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_summary.csv")
meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")

# Merged per-component c2l + NMF
MERGED_CSV = ROOT / "tls_component_features_merged.csv"
if MERGED_CSV.exists():
    w = pd.read_csv(MERGED_CSV)
    print(f"Loaded merged: {len(w)} components")
else:
    # Fallback: use NMF weights only
    w = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
    print(f"Fallback: {len(w)} components (no c2l)")

ct_list = list(basis.index)
N_ECO = 4
eco_list = ["E1","E2","E3","E4"]

eco_names = {"E1":"Lymphocyte TLS","E2":"ILC-enriched TLS",
             "E3":"Myeloid-vascular TLS","E4":"Glial-CD4 TLS"}
eco_short = ["Lymphocyte","ILC-enriched","Myeloid-vascular","Glial-CD4"]
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]
eco_n = [int(summary.loc[i,"n_units"]) for i in range(N_ECO)]
eco_samples = [int(summary.loc[i,"n_samples"]) for i in range(N_ECO)]

if "maturity" not in w.columns and "B" in w.columns:
    ilc_cts = [c for c in w.columns if c.startswith("ILC")]
    mat_cts = [c for c in ["B","Plasma","CD4_T","CD8_T","Dendritic","NK"] if c in w.columns]
    w["ILC_total"] = w[ilc_cts].sum(axis=1) if ilc_cts else np.nan
    w["maturity"] = w[mat_cts].sum(axis=1) if mat_cts else np.nan
w["eco_name"] = w["dominant_ecotype"].map(eco_names)

# ============================================================
# Fig 1: NMF basis heatmap
# ============================================================
basis_z = StandardScaler().fit_transform(basis.values.T).T
vmax = max(abs(basis_z.min()), abs(basis_z.max()))

fig, ax = plt.subplots(figsize=(5, 4.5))
im = ax.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_yticks(range(len(ct_list))); ax.set_yticklabels(ct_list, fontsize=7)
for i, e in enumerate(eco_list):
    ax.text(i, len(ct_list)+0.4, f"{eco_names[e]}\nn={eco_n[i]}", ha="center", va="top",
            fontsize=7, fontweight="bold", color=colors[i])
ax.set_xticks([]); ax.set_xlim(-0.5, N_ECO-0.5); ax.set_ylim(len(ct_list)+1.8, -0.5)
plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02).set_label("z-score", fontsize=7)
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_nmf_basis.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("1/7 fig_nmf_basis")

# ============================================================
# Fig 2: NMF coefficients sorted heatmap
# ============================================================
H = w[eco_list].values
dominant = np.argmax(H, axis=1)
max_w = H.max(axis=1)
order = np.lexsort((max_w, dominant))
H_sorted = H[order]; dom_sorted = dominant[order]
sizes = [(dom_sorted==i).sum() for i in range(N_ECO)]

fig, ax = plt.subplots(figsize=(5, 8))
im = ax.imshow(H_sorted.T, aspect="auto", cmap="YlOrRd", vmin=0, vmax=H_sorted.max()*0.9)
ax.set_yticks(range(N_ECO))
ax.set_yticklabels([f"{eco_short[i]} (n={sizes[i]})" for i in range(N_ECO)], fontsize=8)
for i in range(N_ECO):
    ax.axhline(y=i+0.5, color="black", linewidth=0.5)
boundary = 0
for i in range(N_ECO-1):
    boundary += sizes[i]
    ax.axvline(x=boundary-0.5, color="black", linewidth=0.8, linestyle="--")
for i in range(N_ECO):
    start = sum(sizes[:i]); end = start + sizes[i]
    ax.axvspan(start-0.5, end-0.5, ymin=1.02, ymax=1.06, color=colors[i], clip_on=False, alpha=0.8)
ax.set_xticks([]); ax.set_xlabel(f"{len(w)} TLS components", fontsize=8)
plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02).set_label("NMF weight", fontsize=7)
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_nmf_coefficients.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("2/7 fig_nmf_coefficients")

# ============================================================
# Fig 3: Maturity + ILC boxplots
# ============================================================
rng = np.random.RandomState(42)
def jitter(n, s=0.14): return rng.uniform(-s, s, n)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.8))
for ax, feat, label in [(ax1,"maturity","Maturity score"), (ax2,"ILC_total","ILC total abundance")]:
    data = [w[w["dominant_ecotype"]==e][feat].dropna().values for e in eco_list]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                     medianprops={"color":"black","linewidth":0.8},
                     flierprops={"marker":"none"})
    for i, (patch, vals) in enumerate(zip(bp["boxes"], data)):
        patch.set_facecolor(colors[i]); patch.set_alpha(0.35)
        ax.scatter(np.full(len(vals), i+1) + jitter(len(vals)), vals,
                   s=4, alpha=0.25, color=colors[i], linewidth=0)
    ax.set_xticklabels(eco_short, fontsize=7, rotation=25, ha="right")
    ax.set_ylabel(label, fontsize=8)
fig.tight_layout()
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_ecotype_features.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("3/7 fig_ecotype_features")

# ============================================================
# Fig 4: Maturity vs ILC scatter
# ============================================================
fig, ax = plt.subplots(figsize=(5.5, 4.5))
for i, e in enumerate(eco_list):
    sub = w[w["dominant_ecotype"]==e]
    ax.scatter(sub["maturity"], sub["ILC_total"], s=7, alpha=0.3, color=colors[i],
               label=f"{eco_short[i]} (n={len(sub)})", linewidth=0)
ax.set_xlabel("Maturity score", fontsize=9)
ax.set_ylabel("ILC total abundance", fontsize=9)
ax.legend(fontsize=7, loc="lower right", markerscale=1.2)
fig.tight_layout()
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_maturity_vs_ilc.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("4/7 fig_maturity_vs_ilc")

# ============================================================
# Fig 5: ILC proportions stacked bar
# ============================================================
ilc_pct = w.groupby("dominant_ecotype")[ilc_cts].mean()
total = ilc_pct.sum(axis=1)
ilc_pct_pct = ilc_pct.div(total, axis=0) * 100
ilc_pct_pct = ilc_pct_pct.reindex(eco_list)

fig, ax = plt.subplots(figsize=(5, 4))
bottom = np.zeros(N_ECO)
bar_cols = ["#fbb4ae","#b3cde3","#ccebc5"]
for j, (ct, c) in enumerate(zip(["ILC1","ILC2","ILC3"], bar_cols)):
    vals = ilc_pct_pct[ct].values
    ax.bar(range(N_ECO), vals, bottom=bottom, color=c, alpha=0.8, label=ct,
           edgecolor="black", linewidth=0.3, width=0.55)
    bottom += vals
ax.set_xticks(range(N_ECO)); ax.set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
ax.set_ylabel("% of total ILC", fontsize=9); ax.set_ylim(0, 105)
ax.set_title("ILC subtype proportions (no bias)", fontsize=10, fontweight="bold")
ax.legend(fontsize=8, loc="upper right", ncol=3)
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_ilc_proportions.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("5/7 fig_ilc_proportions")

# ============================================================
# Fig 6: Component & sample counts
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
ax1.bar(range(N_ECO), eco_n, color=colors, alpha=0.7, edgecolor="black", linewidth=0.3)
for i, v in enumerate(eco_n): ax1.text(i, v+1, str(v), ha="center", fontsize=8)
ax1.set_xticks(range(N_ECO)); ax1.set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
ax1.set_ylabel("TLS components", fontsize=9)
ax1.set_title("Components per ecotype", fontsize=10, fontweight="bold")

ax2.bar(range(N_ECO), eco_samples, color=colors, alpha=0.5, edgecolor="black", linewidth=0.3)
for i, v in enumerate(eco_samples): ax2.text(i, v+0.5, str(v), ha="center", fontsize=8)
ax2.set_xticks(range(N_ECO)); ax2.set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
ax2.set_ylabel("Samples", fontsize=9)
ax2.set_title("Sample coverage", fontsize=10, fontweight="bold")
fig.tight_layout()
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_ecotype_counts.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("6/7 fig_ecotype_counts")

# ============================================================
# Fig 7: TLS size boxplot
# ============================================================
# Use metadata for spot counts
meta["dominant_ecotype"] = meta["unit_id"].map(
    dict(zip(w["unit_id"], w["dominant_ecotype"])))
meta_valid = meta.dropna(subset=["dominant_ecotype"])

fig, ax = plt.subplots(figsize=(5, 3.8))
data = [meta_valid[meta_valid["dominant_ecotype"]==e]["n_spots"].values for e in eco_list]
bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                 medianprops={"color":"black","linewidth":0.8}, flierprops={"marker":"none"})
for i, (patch, vals) in enumerate(zip(bp["boxes"], data)):
    patch.set_facecolor(colors[i]); patch.set_alpha(0.35)
    ax.scatter(np.full(len(vals), i+1) + jitter(len(vals)), vals,
               s=4, alpha=0.25, color=colors[i], linewidth=0)
ax.set_xticklabels(eco_short, fontsize=7, rotation=25, ha="right")
ax.set_ylabel("Spots per component", fontsize=8)
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_tls_size.{fmt}", dpi=300 if fmt=="jpg" else None, bbox_inches="tight")
plt.close()
print("7/7 fig_tls_size")

# ====== Summary ======
print("\n=== Ecotype summary (post-fix) ===")
for i, e in enumerate(eco_list):
    top3 = basis[e].sort_values(ascending=False).head(4)
    top3_str = ", ".join([f"{ct}({v:.1f})" for ct, v in top3.items()])
    sub = w[w["dominant_ecotype"]==e]
    print(f"{eco_names[e]}: n={eco_n[i]}, samples={eco_samples[i]}, "
          f"maturity={sub['maturity'].mean():.1f}, ILC={sub['ILC_total'].mean():.2f}")
    print(f"  Basis: {top3_str}")
