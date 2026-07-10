"""
Final integrated figure — classic NMF style:
- Panel A: NMF basis heatmap with top marker genes annotated on columns
- Panel B: GO enrichment bar plots per ecotype
- Panel C: Ecotype features (maturity, ILC, TLS score)
- Panel D: Component & sample counts
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")
OUT = ROOT

# ====== Load data ======
basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_annotated_summary.csv")
go = pd.read_csv(ROOT / "tls_nmf_ecotype_go_enrichment.csv")

ct_list = list(basis.index)
N_ECO = 4
eco_names_list = [f"E{i+1}" for i in range(N_ECO)]

ECO_NAMES = {"E1": "Lymphocyte\nTLS", "E2": "ILC-enriched\nTLS",
             "E3": "Myeloid-vascular\nTLS", "E4": "Glial-CD4\nTLS"}
eco_short = ["Lymphocyte","ILC-enriched","Myeloid-vascular","Glial-CD4"]
palette = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]

# ====== Top genes per ecotype ======
top_genes = {}
for eco in eco_names_list:
    deg = pd.read_csv(ROOT / f"tls_compnmf_rank4_{eco}_top30_up.csv")
    top_genes[eco] = deg["gene"].head(8).tolist()

# ====== GO:BP top terms per ecotype ======
go_bp = go[(go["Gene_set"]=="GO_Biological_Process_2023")].copy()
go_bp["neg_log10_fdr"] = -np.log10(go_bp["Adjusted P-value"].clip(lower=1e-30))

# Top 5 GO:BP per ecotype
go_top = {}
for eco in eco_names_list:
    bp = go_bp[go_bp["ecotype"]==eco].sort_values("Adjusted P-value").head(5)
    go_top[eco] = bp

# ====== Figure ======
fig = plt.figure(figsize=(16, 11))

# ---- Panel A: NMF basis heatmap + gene annotations ----
ax1 = fig.add_axes([0.04, 0.40, 0.42, 0.56])
basis_z = StandardScaler().fit_transform(basis.values.T).T
vmax = max(abs(basis_z.min()), abs(basis_z.max()))
im1 = ax1.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

# Cell type labels
ax1.set_yticks(range(len(ct_list)))
ax1.set_yticklabels(ct_list, fontsize=6.5)

# Ecotype labels on top
eco_n = [summary.loc[i, "n_units"] for i in range(N_ECO)]
for i in range(N_ECO):
    ax1.text(i, -0.8, ECO_NAMES[eco_names_list[i]], ha="center", va="top", fontsize=7,
             fontweight="bold", color=palette[i])
    ax1.text(i, -0.1, f"n={eco_n[i]}", ha="center", va="top", fontsize=6, color="grey")

# Annotate top genes below each ecotype column
for i, eco in enumerate(eco_names_list):
    genes = top_genes[eco][:5]
    for j, g in enumerate(genes):
        y_pos = len(ct_list) + 0.5 + j * 0.55
        ax1.text(i, y_pos, g, ha="center", va="center", fontsize=4.8, style="italic",
                 color=palette[i], bbox={"boxstyle":"round,pad=0.15","fc":"white","ec":palette[i],"lw":0.3,"alpha":0.8})

ax1.set_xticks([])
ax1.set_xlim(-0.5, N_ECO - 0.5)
ax1.set_ylim(len(ct_list) + 3.2, -1.5)
ax1.set_title("a  NMF basis (cell type weights, z-score) + top marker genes", fontsize=9, fontweight="bold", loc="left")
plt.colorbar(im1, ax=ax1, shrink=0.6, pad=0.02).set_label("z-score", fontsize=6)

# ---- Panel B: GO enrichment bar plots ----
ax2 = fig.add_axes([0.52, 0.40, 0.46, 0.56])

# Collect all terms
all_terms = []
all_fdr = []
all_eco = []
for eco in eco_names_list:
    bp = go_top[eco]
    for _, row in bp.iterrows():
        term = row["Term"]
        # Clean GO prefix and truncate
        term = term.replace(" (GO:", "|").split("|")[0].strip()
        if len(term) > 55:
            term = term[:52] + "..."
        all_terms.append(term)
        all_fdr.append(row["neg_log10_fdr"])
        all_eco.append(eco)

y_positions = list(range(len(all_terms)))[::-1]
colors = [palette[eco_names_list.index(e)] for e in all_eco]
ax2.barh(y_positions, all_fdr, color=colors, alpha=0.7, edgecolor="black", linewidth=0.3, height=0.7)
ax2.set_yticks(y_positions)
ax2.set_yticklabels(all_terms, fontsize=5.5)
ax2.set_xlabel("-log10(FDR)", fontsize=7)
ax2.set_title("b  GO Biological Process enrichment (top 5 per ecotype)", fontsize=9, fontweight="bold", loc="left")

# Add ecotype group labels
eco_boundaries = {}
for i, eco in enumerate(eco_names_list):
    eco_indices = [j for j, e in enumerate(all_eco) if e == eco]
    if eco_indices:
        eco_boundaries[eco] = (min(eco_indices), max(eco_indices))
# Draw group brackets on right side
for i, eco in enumerate(eco_names_list):
    if eco in eco_boundaries:
        mn, mx = eco_boundaries[eco]
        mid = (y_positions[mn] + y_positions[mx]) / 2
        ax2.text(ax2.get_xlim()[1] * 0.98, mid, ECO_NAMES[eco].replace("\n"," "),
                 fontsize=6.5, fontweight="bold", color=palette[i], va="center", ha="right")

# ---- Panel C: Ecotype features ----
ax3 = fig.add_axes([0.04, 0.05, 0.42, 0.28])
feat_mat = summary[["maturity_score","ILC_total","mean_tls_score"]].values
feat_names = ["Maturity","ILC total\nabundance","TLS score"]
feat_z = (feat_mat - feat_mat.mean(axis=0)) / (feat_mat.std(axis=0) + 1e-8)
im3 = ax3.imshow(feat_z.T, aspect="auto", cmap="YlOrRd", vmin=-1.2, vmax=1.5)
ax3.set_yticks(range(len(feat_names)))
ax3.set_yticklabels(feat_names, fontsize=6.5)
ax3.set_xticks(range(N_ECO))
ax3.set_xticklabels(eco_short, fontsize=6.5, rotation=20, ha="right")
ax3.set_title("c  Ecotype features (z-score)", fontsize=9, fontweight="bold", loc="left")
# Annotate values
for i in range(N_ECO):
    for j in range(len(feat_names)):
        ax3.text(i, j, f"{feat_mat[i,j]:.1f}", ha="center", va="center", fontsize=6,
                 color="white" if abs(feat_z[i,j]) > 0.7 else "black")
plt.colorbar(im3, ax=ax3, shrink=0.7).set_label("z-score", fontsize=6)

# ---- Panel D1: Component counts ----
ax4 = fig.add_axes([0.52, 0.18, 0.22, 0.15])
eco_counts = [summary.loc[i, "n_units"] for i in range(N_ECO)]
ax4.bar(range(N_ECO), eco_counts, color=palette, alpha=0.7, edgecolor="black", linewidth=0.3)
ax4.set_xticks(range(N_ECO))
ax4.set_xticklabels(eco_short, fontsize=6, rotation=20, ha="right")
ax4.set_ylabel("Components", fontsize=7)
ax4.set_title("d  Components per ecotype", fontsize=9, fontweight="bold", loc="left")
for i, v in enumerate(eco_counts):
    ax4.text(i, v + max(eco_counts)*0.02, str(v), ha="center", fontsize=7)

# ---- Panel D2: Sample counts ----
ax5 = fig.add_axes([0.78, 0.18, 0.20, 0.15])
eco_samples = [summary.loc[i, "n_samples"] for i in range(N_ECO)]
ax5.bar(range(N_ECO), eco_samples, color=palette, alpha=0.5, edgecolor="black", linewidth=0.3)
ax5.set_xticks(range(N_ECO))
ax5.set_xticklabels(eco_short, fontsize=6, rotation=20, ha="right")
ax5.set_ylabel("Samples", fontsize=7)
ax5.set_title("e  Sample coverage", fontsize=9, fontweight="bold", loc="left")
for i, v in enumerate(eco_samples):
    ax5.text(i, v + max(eco_samples)*0.02, str(v), ha="center", fontsize=7)

# ---- Panel E: ILC proportions (stacked bar, confirm no bias) ----
ax6 = fig.add_axes([0.52, 0.05, 0.46, 0.09])
ilc_pct = summary[["ILC1_frac","ILC2_frac","ILC3_frac"]].values * 100
bottom = np.zeros(N_ECO)
ilc_colors = ["#fbb4ae","#b3cde3","#ccebc5"]
for j, (ct, c) in enumerate(zip(["ILC1","ILC2","ILC3"], ilc_colors)):
    ax6.bar(range(N_ECO), ilc_pct[:, j], bottom=bottom, color=c, alpha=0.8, label=ct,
            edgecolor="black", linewidth=0.2, width=0.6)
    bottom += ilc_pct[:, j]
ax6.set_xticks(range(N_ECO))
ax6.set_xticklabels(eco_short, fontsize=6, rotation=20, ha="right")
ax6.set_ylabel("% ILC", fontsize=6)
ax6.set_ylim(0, 105)
ax6.set_title("f  ILC subtype proportion (no bias: ~33:33:34)", fontsize=8, fontweight="bold", loc="left")
ax6.legend(fontsize=5.5, loc="upper right", ncol=3)

fig.savefig(OUT / "fig_tls_nmf_ecotypes_final.jpg", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "fig_tls_nmf_ecotypes_final.pdf", bbox_inches="tight")
plt.close()
print(f"Figure: {OUT}/fig_tls_nmf_ecotypes_final.jpg")

# ====== Summary ======
print("\nEcotype summary:")
for i, eco in enumerate(eco_names_list):
    name = ECO_NAMES[eco].replace("\n"," ")
    n = summary.loc[i, "n_units"]
    ns = summary.loc[i, "n_samples"]
    mat = summary.loc[i, "maturity_score"]
    ilc = summary.loc[i, "ILC_total"]
    genes = ", ".join(top_genes[eco][:5])
    go_top_term = go_top[eco].iloc[0]["Term"] if len(go_top[eco]) > 0 else "N/A"
    print(f"  {name}: n={n}, samples={ns}, maturity={mat:.1f}, ILC={ilc:.2f}")
    print(f"    Genes: {genes}")
    print(f"    GO: {go_top_term}")
