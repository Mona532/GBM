"""
All standalone figures — one figure per scientific question.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":8,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")

# ====== Shared data ======
basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_annotated_summary.csv")
go = pd.read_csv(ROOT / "tls_nmf_ecotype_go_enrichment.csv")
go_bp = go[(go["Gene_set"]=="GO_Biological_Process_2023")].copy()
go_bp["neg_log10_fdr"] = -np.log10(go_bp["Adjusted P-value"].clip(lower=1e-30))

N_ECO = 4
eco_labels = ["Lymphocyte\nTLS", "ILC-enriched\nTLS", "Myeloid-vascular\nTLS", "Glial-CD4\nTLS"]
eco_short = ["Lymphocyte", "ILC-enriched", "Myeloid-vascular", "Glial-CD4"]
colors = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]

basis_z = StandardScaler().fit_transform(basis.values.T).T
vmax = max(abs(basis_z.min()), abs(basis_z.max()))
ct_list = list(basis.index)

eco_n = [int(summary.loc[i, "n_units"]) for i in range(N_ECO)]
eco_samples = [int(summary.loc[i, "n_samples"]) for i in range(N_ECO)]
maturity = summary["maturity_score"].values
ilc_total = summary["ILC_total"].values
tls_score = summary["mean_tls_score"].values
ilc_fracs = summary[["ILC1_frac","ILC2_frac","ILC3_frac"]].values

# ============================================================
# Fig 1: NMF basis heatmap
# ============================================================
fig, ax = plt.subplots(figsize=(5, 4.5))
im = ax.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_yticks(range(len(ct_list))); ax.set_yticklabels(ct_list, fontsize=7)
for i in range(N_ECO):
    ax.text(i, len(ct_list)+0.4, eco_labels[i], ha="center", va="top", fontsize=7, fontweight="bold", color=colors[i])
    ax.text(i, len(ct_list)+2.1, f"n={eco_n[i]}", ha="center", va="top", fontsize=6.5, color="grey")
ax.set_xticks([]); ax.set_xlim(-0.5, N_ECO-0.5); ax.set_ylim(len(ct_list)+2.8, -0.5)
cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02); cbar.set_label("z-score", fontsize=7)
fig.savefig(ROOT / "fig_nmf_basis.jpg", dpi=300, bbox_inches="tight")
fig.savefig(ROOT / "fig_nmf_basis.pdf", bbox_inches="tight")
plt.close()
print("1/7 fig_nmf_basis.jpg")

# ============================================================
# Fig 2: GO enrichment bar plot
# ============================================================
fig, ax = plt.subplots(figsize=(8, 5))
all_terms, all_fdr, all_colors = [], [], []
for i, eco_idx in enumerate(range(N_ECO)):
    eco = f"E{eco_idx+1}"
    bp = go_bp[go_bp["ecotype"]==eco].sort_values("Adjusted P-value").head(5)
    for _, row in bp.iterrows():
        term = row["Term"]
        term = term.replace(" (GO:", "|").split("|")[0].strip()
        if len(term) > 60: term = term[:57]+"..."
        all_terms.append(term)
        all_fdr.append(row["neg_log10_fdr"])
        all_colors.append(colors[i])

y = list(range(len(all_terms)))[::-1]
ax.barh(y, all_fdr, color=all_colors, alpha=0.75, edgecolor="black", linewidth=0.3, height=0.7)
ax.set_yticks(y); ax.set_yticklabels(all_terms, fontsize=6.5)
ax.set_xlabel("-log10(FDR)", fontsize=8)
ax.set_title("GO Biological Process (top 5 per ecotype)", fontsize=10, fontweight="bold")
# Legend
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=colors[i], label=eco_short[i]) for i in range(N_ECO)],
          fontsize=7, loc="lower right")
fig.savefig(ROOT / "fig_go_enrichment.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("2/7 fig_go_enrichment.jpg")

# ============================================================
# Fig 3: Ecotype features (maturity, ILC, TLS score)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
feats = [("Maturity score", maturity), ("ILC total abundance", ilc_total), ("TLS score", tls_score)]
for ax, (name, vals) in zip(axes, feats):
    ax.bar(range(N_ECO), vals, color=colors, alpha=0.7, edgecolor="black", linewidth=0.3)
    ax.set_xticks(range(N_ECO)); ax.set_xticklabels(eco_short, fontsize=7, rotation=20, ha="right")
    ax.set_ylabel(name, fontsize=8)
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals)*0.02, f"{v:.1f}", ha="center", fontsize=7)
fig.tight_layout()
fig.savefig(ROOT / "fig_ecotype_features.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("3/7 fig_ecotype_features.jpg")

# ============================================================
# Fig 4: ILC proportions stacked bar
# ============================================================
fig, ax = plt.subplots(figsize=(5, 4))
ilc_pct = ilc_fracs * 100
bottom = np.zeros(N_ECO)
ilc_colors_bar = ["#fbb4ae","#b3cde3","#ccebc5"]
for j, (ct, c) in enumerate(zip(["ILC1","ILC2","ILC3"], ilc_colors_bar)):
    ax.bar(range(N_ECO), ilc_pct[:, j], bottom=bottom, color=c, alpha=0.8, label=ct,
           edgecolor="black", linewidth=0.3, width=0.55)
    bottom += ilc_pct[:, j]
ax.set_xticks(range(N_ECO)); ax.set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
ax.set_ylabel("% of total ILC", fontsize=9)
ax.set_ylim(0, 105)
ax.set_title("ILC subtype proportions (no bias)", fontsize=10, fontweight="bold")
ax.legend(fontsize=8, loc="upper right", ncol=3)
fig.savefig(ROOT / "fig_ilc_proportions.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("4/7 fig_ilc_proportions.jpg")

# ============================================================
# Fig 5: TLS counts
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
axes[0].bar(range(N_ECO), eco_n, color=colors, alpha=0.7, edgecolor="black", linewidth=0.3)
axes[0].set_xticks(range(N_ECO)); axes[0].set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
axes[0].set_ylabel("TLS components", fontsize=9)
axes[0].set_title("Components per ecotype", fontsize=10, fontweight="bold")
for i, v in enumerate(eco_n):
    axes[0].text(i, v+1, str(v), ha="center", fontsize=8)

axes[1].bar(range(N_ECO), eco_samples, color=colors, alpha=0.5, edgecolor="black", linewidth=0.3)
axes[1].set_xticks(range(N_ECO)); axes[1].set_xticklabels(eco_short, fontsize=8, rotation=20, ha="right")
axes[1].set_ylabel("Samples", fontsize=9)
axes[1].set_title("Sample coverage", fontsize=10, fontweight="bold")
for i, v in enumerate(eco_samples):
    axes[1].text(i, v+0.5, str(v), ha="center", fontsize=8)

fig.tight_layout()
fig.savefig(ROOT / "fig_ecotype_counts.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("5/7 fig_ecotype_counts.jpg")

# ============================================================
# Fig 6: BANKSY niche composition
# ============================================================
bn = pd.read_csv(ROOT / "banksy_niche_tls_global.csv")
bn_ct = bn[ct_list].values
bn_tls = bn["tls_log2fc"].values
K_B = len(bn)
vmax_b = max(abs(bn_ct.min()), abs(bn_ct.max()), 0.5)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw={"width_ratios":[2,1]})
im = axes[0].imshow(bn_ct.T, aspect="auto", cmap="RdBu_r", vmin=-vmax_b, vmax=vmax_b)
axes[0].set_yticks(range(len(ct_list))); axes[0].set_yticklabels(ct_list, fontsize=7)
axes[0].set_xticks(range(K_B))
axes[0].set_xticklabels([f"N{k}\nn={bn.loc[k,'n_spots']:.0f}" for k in range(K_B)], fontsize=7)
axes[0].set_title("BANKSY niche composition", fontsize=10, fontweight="bold")
plt.colorbar(im, ax=axes[0], shrink=0.7).set_label("log2 enrichment", fontsize=7)

tlsc = ["#d7191c" if v>0 else "#2b83ba" for v in bn_tls]
axes[1].bar(range(K_B), bn_tls, color=tlsc, alpha=0.7, edgecolor="black", linewidth=0.3)
axes[1].axhline(0, color="black", linewidth=0.5, linestyle="--")
axes[1].set_xticks(range(K_B)); axes[1].set_xticklabels([f"N{k}" for k in range(K_B)], fontsize=8)
axes[1].set_ylabel("TLS enrichment (log2 FC)", fontsize=8)
axes[1].set_title("TLS enrichment", fontsize=10, fontweight="bold")
for i, v in enumerate(bn_tls):
    axes[1].text(i, v+0.05*np.sign(v), f"{v:+.2f}", ha="center", fontsize=8)

fig.tight_layout()
fig.savefig(ROOT / "fig_banksy_niche.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("6/7 fig_banksy_niche.jpg")

# ============================================================
# Fig 7: Maturity vs ILC scatter
# ============================================================
comp_w = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
comp_d = pd.read_csv(ROOT / "tls_components_denovo.csv")
comp = comp_w.merge(comp_d[["sample","component_id","B_mean","Plasma_mean","CD4_T_mean",
    "CD8_T_mean","Dendritic_mean","NK_mean","ILC1_mean","ILC2_mean","ILC3_mean"]],
    on=["sample","component_id"], how="left")
comp["ILC_total"] = comp["ILC1_mean"]+comp["ILC2_mean"]+comp["ILC3_mean"]
comp["maturity"] = (comp["B_mean"]+comp["Plasma_mean"]+comp["CD4_T_mean"]+
                     comp["CD8_T_mean"]+comp["Dendritic_mean"]+comp["NK_mean"])
comp = comp.dropna(subset=["ILC_total","maturity"])

fig, ax = plt.subplots(figsize=(6, 5))
for i in range(N_ECO):
    eco = f"E{i+1}"
    sub = comp[comp["dominant_ecotype"]==eco]
    ax.scatter(sub["maturity"], sub["ILC_total"], s=6, alpha=0.3, color=colors[i],
               label=f"{eco_short[i]} (n={len(sub)})")
ax.set_xlabel("Maturity score", fontsize=9)
ax.set_ylabel("ILC total abundance", fontsize=9)
ax.set_title("Maturity vs ILC per TLS component", fontsize=10, fontweight="bold")
ax.legend(fontsize=7, loc="lower right", markerscale=1.5)
fig.savefig(ROOT / "fig_maturity_vs_ilc.jpg", dpi=300, bbox_inches="tight")
plt.close()
print("7/7 fig_maturity_vs_ilc.jpg")

print("\nAll figures saved to results/")
