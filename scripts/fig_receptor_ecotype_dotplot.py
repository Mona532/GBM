"""Receptor dotplot: 5 categories x 4 NMF ecotypes — both pseudobulk and spot versions"""
import pandas as pd, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import Normalize

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],
    "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,
    "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.6,"legend.frameon":False})

ROOT = "E:/GBM/results"
rx = pd.read_excel("E:/GBM/ti tianran2.xlsx")

cat_names = ["Glutamate","GABA/Gly","Cholinergic","DA/NE","Serotonin"]
cat_colors = ["#E64A19","#1B5E20","#0D47A1","#6A1B9A","#BF360C"]
eco_names = {"E1":"Lymphocyte","E2":"ILC-enriched","E3":"Myeloid-vascular","E4":"Glial-CD4"}
eco_list = ["E1","E2","E3","E4"]

cat_map = {}
for ci, col in enumerate(rx.columns):
    for g in rx[col].dropna():
        cat_map[g] = cat_names[ci]

datasets = [
    ("receptor_pseudobulk_ecotype.csv", "pseudobulk"),
    ("receptor_ecotype_detection.csv", "spot"),
]

for csv_file, suffix in datasets:
    df = pd.read_csv(f"{ROOT}/{csv_file}")
    df["category"] = df["gene"].map(cat_map)
    df = df.dropna(subset=["category"])

    col_groups = []
    max_n = 0
    for cat in cat_names:
        for eco in eco_list:
            sub = df[(df["category"]==cat)&(df["ecotype"]==eco)].sort_values("detect_rate", ascending=False)
            col_groups.append((cat, eco, sub))
            max_n = max(max_n, len(sub))

    n_cols = len(col_groups)
    fig, ax = plt.subplots(figsize=(18, max_n*0.20+1.2))

    vmax = df["mean_expr"].quantile(0.95)
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = mpl.colormaps["YlOrRd"]

    for ci, (cat, eco, sub) in enumerate(col_groups):
        first_in_cat = (ci % 4 == 0)  # only first ecotype per category
        for i, (_, r) in enumerate(sub.iterrows()):
            s = r["detect_rate"]*140+6
            c = cmap(norm(r["mean_expr"]))
            ax.scatter(ci, i, s=s, c=[c], alpha=0.85, edgecolors="white", linewidths=0.3)
            if first_in_cat:
                ax.text(ci+0.42, i, r["gene"], fontsize=5, va="center")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([eco_names[e] for _, e, _ in col_groups], fontsize=5.5, rotation=30, ha="right")
    for ci in range(5):
        mid = ci*4+1.5
        ax.text(mid, -1.5, cat_names[ci], fontsize=7, fontweight="bold", color=cat_colors[ci], ha="center", va="bottom")
        if ci > 0:
            ax.axvline(x=ci*4-0.5, color="lightgrey", linewidth=0.5, linestyle="--")

    ax.set_yticks([]); ax.set_xlim(-0.8, n_cols-0.2); ax.invert_yaxis()

    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.25, aspect=16, pad=0.02)
    cbar.set_label("Mean CPM", fontsize=6); cbar.ax.tick_params(labelsize=5)

    for pct, label, xoff in [(25,"25%",-0.65),(50,"50%",-0.50),(75,"75%",-0.30)]:
        s = pct/100*140+6
        ax.scatter(xoff, -1.0, s=s, c="grey", alpha=0.35, edgecolors="black", linewidths=0.3, clip_on=False)
        ax.text(xoff, -0.3, label, fontsize=5, ha="center")

    fig.tight_layout()
    for fmt in ["jpg","svg","pdf"]:
        fig.savefig(f"{ROOT}/fig_receptor_ecotype_dotplot_{suffix}.{fmt}",
                    dpi=600 if fmt=="jpg" else None, bbox_inches="tight")
    plt.close()
    print(f"Done: fig_receptor_ecotype_dotplot_{suffix}")

print("All done.")
