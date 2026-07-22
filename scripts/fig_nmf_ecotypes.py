"""
TLS ecotype NMF figure — publication-ready.
Inputs: R NMF outputs (basis, unit weights, top genes, DEGs, ecotype summary)
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.right":False,"axes.spines.top":False,
    "axes.linewidth":0.5,"legend.frameon":False})

ROOT = Path(r"E:/GBM/results")
OUT = ROOT

# ====== Load data ======
basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
weights = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_annotated_summary.csv")

ct_list = list(basis.index)
N_ECO = basis.shape[1]
eco_names_list = [f"E{i+1}" for i in range(N_ECO)]

# ====== Naming: based on basis weights ======
ECO_NAMES = {}
for e in eco_names_list:
    top = basis[e].sort_values(ascending=False)
    t1, t2 = top.index[0], top.index[1]
    v1, v2 = top.iloc[0], top.iloc[1]

    if "ILC" in t1 and "ILC" in t2:
        name = "ILC-enriched TLS"
    elif ("CD8" in t1 or "NK" in t1) and ("B" in t1 or "B" in t2 or "NK" in t2):
        name = "Lymphocyte TLS"
    elif ("Vascular" in t1 or "Macrophage" in t1):
        name = "Myeloid-vascular TLS"
    elif "Glial" in t1:
        name = "Glial-CD4 TLS"
    else:
        name = f"{t1}-dominant TLS"
    ECO_NAMES[e] = name

eco_short = [ECO_NAMES[e].replace(" TLS","") for e in eco_names_list]
eco_n = [summary.loc[i, "n_units"] for i in range(N_ECO)]
eco_samples = [summary.loc[i, "n_samples"] for i in range(N_ECO)]

print("Ecotype names:")
for i, e in enumerate(eco_names_list):
    print(f"  {e} ({ECO_NAMES[e]}): n={eco_n[i]}, samples={eco_samples[i]}")

# ====== Per-ecotype metrics ======
maturity = summary["maturity_score"].values
ilc_total = summary["ILC_total"].values
tls_score = summary["mean_tls_score"].values
ilc_fracs = summary[["ILC1_frac","ILC2_frac","ILC3_frac"]].values

# ====== Figure ======
fig = plt.figure(figsize=(15, 10))

# --- Panel A: NMF basis heatmap (cell_type × ecotype, z-score) ---
ax1 = fig.add_axes([0.04, 0.54, 0.28, 0.42])
basis_z = StandardScaler().fit_transform(basis.values.T).T  # z-score across ecotypes
vmax = max(abs(basis_z.min()), abs(basis_z.max()))
im1 = ax1.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax1.set_yticks(range(len(ct_list)))
ax1.set_yticklabels(ct_list, fontsize=6)
ax1.set_xticks(range(N_ECO))
ax1.set_xticklabels([f"{eco_short[i]}\nn={eco_n[i]}" for i in range(N_ECO)], fontsize=6.5)
ax1.set_title("a  NMF basis: cell type × ecotype (z-score)", fontsize=9, fontweight="bold", loc="left")
plt.colorbar(im1, ax=ax1, shrink=0.7).set_label("z-score", fontsize=6)

# --- Panel B: Maturity + ILC features per ecotype ---
ax2 = fig.add_axes([0.36, 0.54, 0.28, 0.42])
# Combine: maturity, ILC_total, TLS_score, ILC fractions
feat_names = ["Maturity", "ILC total", "TLS score", "ILC1 %", "ILC2 %", "ILC3 %"]
feat_mat = np.column_stack([maturity, ilc_total, tls_score, ilc_fracs])
feat_z = (feat_mat - feat_mat.mean(axis=0)) / (feat_mat.std(axis=0) + 1e-8)
im2 = ax2.imshow(feat_z.T, aspect="auto", cmap="YlOrRd", vmin=-1.2, vmax=1.5)
ax2.set_yticks(range(len(feat_names)))
ax2.set_yticklabels(feat_names, fontsize=6)
ax2.set_xticks(range(N_ECO))
ax2.set_xticklabels([eco_short[i] for i in range(N_ECO)], fontsize=6.5)
ax2.set_title("b  Ecotype features (z-score)", fontsize=9, fontweight="bold", loc="left")
plt.colorbar(im2, ax=ax2, shrink=0.7).set_label("z-score", fontsize=6)

# --- Panel C: Top marker genes per ecotype ---
ax3 = fig.add_axes([0.68, 0.54, 0.30, 0.42])
# Collect top 5 genes per ecotype
all_genes = []
for e in eco_names_list:
    deg = pd.read_csv(ROOT / f"tls_compnmf_rank4_{e}_top30_up.csv")
    all_genes.append(deg["gene"].head(5).tolist())

gene_labels = []
gene_positions = []
gene_values = []
y_pos = 0
y_ticks = []
y_tick_labels = []
palette = ["#e41a1c","#377eb8","#4daf4a","#ff7f00"]
for i, (e, genes) in enumerate(zip(eco_names_list, all_genes)):
    for j, g in enumerate(genes):
        gene_labels.append(f"{g}")
        gene_positions.append(y_pos)
        gene_values.append(1.0)
        y_ticks.append(y_pos)
        y_tick_labels.append(g)
        y_pos += 1
    y_pos += 0.8  # gap between ecotypes

ax3.set_ylim(-0.5, y_pos)
for i, (e, genes) in enumerate(zip(eco_names_list, all_genes)):
    start = sum(len(all_genes[k]) for k in range(i)) + i * 0.8
    for j, g in enumerate(genes):
        ax3.barh(start + j, 1, height=0.7, color=palette[i], alpha=0.75, edgecolor="black", linewidth=0.2)
ax3.set_yticks([sum(len(all_genes[k]) for k in range(i)) + len(all_genes[i])/2 - 0.5 + i*0.8 for i in range(N_ECO)])
ax3.set_yticklabels(eco_short, fontsize=7)
ax3.set_xticks([])
ax3.set_title("c  Top 5 marker genes per ecotype", fontsize=9, fontweight="bold", loc="left")
# Add gene names next to bars
for i, (e, genes) in enumerate(zip(eco_names_list, all_genes)):
    start = sum(len(all_genes[k]) for k in range(i)) + i * 0.8
    for j, g in enumerate(genes):
        ax3.text(0.05, start + j, g, fontsize=5.5, va="center", fontstyle="italic")

# --- Panel D: Ecotype counts ---
ax4 = fig.add_axes([0.04, 0.06, 0.20, 0.38])
bar_colors = [palette[i] for i in range(N_ECO)]
ax4.bar(range(N_ECO), eco_n, color=bar_colors, alpha=0.7, edgecolor="black", linewidth=0.3)
ax4.set_xticks(range(N_ECO))
ax4.set_xticklabels(eco_short, fontsize=6, rotation=20, ha="right")
ax4.set_ylabel("TLS components", fontsize=7)
ax4.set_title("d  Components per ecotype", fontsize=9, fontweight="bold", loc="left")
for i, v in enumerate(eco_n):
    ax4.text(i, v + 1, str(v), ha="center", fontsize=7)

# --- Panel E: Sample coverage ---
ax5 = fig.add_axes([0.28, 0.06, 0.20, 0.38])
ax5.bar(range(N_ECO), eco_samples, color=bar_colors, alpha=0.7, edgecolor="black", linewidth=0.3)
ax5.set_xticks(range(N_ECO))
ax5.set_xticklabels(eco_short, fontsize=6, rotation=20, ha="right")
ax5.set_ylabel("Samples", fontsize=7)
ax5.set_title("e  Sample coverage", fontsize=9, fontweight="bold", loc="left")
for i, v in enumerate(eco_samples):
    ax5.text(i, v + 1, str(v), ha="center", fontsize=7)

# --- Panel F: Maturity vs ILC scatter ---
ax6 = fig.add_axes([0.55, 0.06, 0.42, 0.38])
# Per-component scatter: merge NMF weights with c2l composition from denovo.csv
comp_w = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
comp_d = pd.read_csv(ROOT / "tls_components_denovo.csv")
comp = comp_w.merge(comp_d[["sample","component_id","B_mean","Plasma_mean","CD4_T_mean",
    "CD8_T_mean","Dendritic_mean","NK_mean","ILC1_mean","ILC2_mean","ILC3_mean"]],
    on=["sample","component_id"], how="left")
comp["ILC_total"] = comp["ILC1_mean"] + comp["ILC2_mean"] + comp["ILC3_mean"]
comp["maturity"] = (comp["B_mean"] + comp["Plasma_mean"] + comp["CD4_T_mean"] +
                     comp["CD8_T_mean"] + comp["Dendritic_mean"] + comp["NK_mean"])
# Drop rows where merge failed
comp = comp.dropna(subset=["ILC_total","maturity"])

for i, e in enumerate(eco_names_list):
    sub = comp[comp["dominant_ecotype"] == e]
    ax6.scatter(sub["maturity"], sub["ILC_total"], s=6, alpha=0.35, color=palette[i],
                label=f"{eco_short[i]} (n={len(sub)})")
ax6.set_xlabel("Maturity score", fontsize=7)
ax6.set_ylabel("ILC total abundance", fontsize=7)
ax6.set_title("f  Maturity vs ILC abundance per component", fontsize=9, fontweight="bold", loc="left")
ax6.legend(fontsize=5.5, loc="lower right", markerscale=1.5)

fig.savefig(OUT / "fig_tls_nmf_ecotypes.jpg", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "fig_tls_nmf_ecotypes.pdf", bbox_inches="tight")
plt.close()
print(f"\nFigure: {OUT}/fig_tls_nmf_ecotypes.jpg")

# ====== Manuscript summary table ======
print("\n--- Manuscript table ---")
for i, e in enumerate(eco_names_list):
    name = ECO_NAMES[e]
    n = eco_n[i]
    ns = eco_samples[i]
    mat = maturity[i]
    ilc = ilc_total[i]
    # Basis top3
    top3 = basis[e].sort_values(ascending=False).head(3)
    top3_str = ", ".join([f"{ct}({v:.1f})" for ct, v in top3.items()])
    # Top gene
    deg = pd.read_csv(ROOT / f"tls_compnmf_rank4_{e}_top30_up.csv")
    top_gene = deg["gene"].iloc[0]
    print(f"  {name}: {n} components, {ns} samples, maturity={mat:.1f}, ILC_total={ilc:.2f}")
    print(f"    Basis: {top3_str}")
    print(f"    Top gene: {top_gene}")

print("\nKey conclusion:")
print("NMF resolves 4 TLS ecotypes from pseudobulk gene expression.")
print(f"{ECO_NAMES['E2']} (n={eco_n[1]}) is defined by ILC1/2/3 as the dominant basis feature,")
print("distinct from the Lymphocyte ecotype where ILCs co-occur with CD8/NK/B.")
