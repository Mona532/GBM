"""
Nature-style boxplots + jittered points for ecotype features.
Replaces bar plots with distribution-level visualization.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.6,
    "legend.frameon": False,
})

ROOT = Path(r"E:/GBM/results")

# ====== Load per-component data ======
df = pd.read_csv(ROOT / "tls_component_features_per_component.csv")
df["ecotype_name"] = df["dominant_ecotype"].map({
    "E1": "Lymphocyte", "E2": "ILC-enriched",
    "E3": "Myeloid-vascular", "E4": "Glial-CD4"
})
df = df.dropna(subset=["ILC_total","maturity_score"])

eco_order = ["Lymphocyte","ILC-enriched","Myeloid-vascular","Glial-CD4"]
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]

def jitter(a, spread=0.15):
    return np.random.RandomState(42).uniform(-spread, spread, len(a))

# ============================================================
# Fig: Maturity + ILC total boxplots
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.8))

for ax, feat, label in [(ax1, "maturity_score", "Maturity score"),
                          (ax2, "ILC_total", "ILC total abundance")]:
    data = [df[df["ecotype_name"]==e][feat].values for e in eco_order]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                     medianprops={"color":"black","linewidth":0.8},
                     flierprops={"marker":"none"},
                     whiskerprops={"linewidth":0.6},
                     capprops={"linewidth":0.6})
    for i, (patch, vals) in enumerate(zip(bp["boxes"], data)):
        patch.set_facecolor(colors[i]); patch.set_alpha(0.35)
        ax.scatter(np.full(len(vals), i+1) + jitter(vals, 0.14), vals,
                   s=4, alpha=0.25, color=colors[i], linewidth=0)
    ax.set_xticklabels(eco_order, fontsize=7, rotation=25, ha="right")
    ax.set_ylabel(label, fontsize=8)

fig.tight_layout()
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_ecotype_features.{fmt}", dpi=300 if fmt=="jpg" else None,
                bbox_inches="tight")
plt.close()
print(f"1/3 fig_ecotype_features (svg/pdf/jpg)")

# ============================================================
# Fig: Maturity vs ILC scatter (already fine, just polish)
# ============================================================
fig, ax = plt.subplots(figsize=(5.5, 4.5))
for i, e in enumerate(eco_order):
    sub = df[df["ecotype_name"]==e]
    ax.scatter(sub["maturity_score"], sub["ILC_total"], s=7, alpha=0.3, color=colors[i],
               label=f"{e} (n={len(sub)})", linewidth=0)
ax.set_xlabel("Maturity score", fontsize=9)
ax.set_ylabel("ILC total abundance", fontsize=9)
ax.legend(fontsize=7, loc="lower right", markerscale=1.2, handletextpad=0.5)
fig.tight_layout()
for fmt in ["jpg","svg","pdf"]:
    fig.savefig(ROOT / f"fig_maturity_vs_ilc.{fmt}", dpi=300 if fmt=="jpg" else None,
                bbox_inches="tight")
plt.close()
print(f"2/3 fig_maturity_vs_ilc (svg/pdf/jpg)")

# ============================================================
# Fig: TLS size boxplot per ecotype
# ============================================================
meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")
df2 = df.merge(meta[["unit_id","n_spots"]], on="unit_id", how="left")

fig, ax = plt.subplots(figsize=(5, 3.8))
data = [df2[df2["ecotype_name"]==e]["n_spots"].dropna().values for e in eco_order]
if all(len(d) > 0 for d in data):
    bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                     medianprops={"color":"black","linewidth":0.8},
                     flierprops={"marker":"none"})
    for i, (patch, vals) in enumerate(zip(bp["boxes"], data)):
        patch.set_facecolor(colors[i]); patch.set_alpha(0.35)
        ax.scatter(np.full(len(vals), i+1) + jitter(vals, 0.14), vals,
                   s=4, alpha=0.25, color=colors[i], linewidth=0)
    ax.set_xticklabels(eco_order, fontsize=7, rotation=25, ha="right")
    ax.set_ylabel("Spots per component", fontsize=8)
    ax.set_ylabel("Spots per component", fontsize=8)
    fig.tight_layout()
    for fmt in ["jpg","svg","pdf"]:
        fig.savefig(ROOT / f"fig_tls_size.{fmt}", dpi=300 if fmt=="jpg" else None,
                    bbox_inches="tight")
    plt.close()
    print(f"3/3 fig_tls_size (svg/pdf/jpg)")

print("\nDone — all boxplots saved.")
