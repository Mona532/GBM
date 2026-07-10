"""Nature-style receptor dotplot: 5 columns, color=expression, size=detection"""
import pandas as pd, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 7, "axes.spines.right": False, "axes.spines.top": False,
    "axes.linewidth": 0.6, "legend.frameon": False,
})

ROOT = "E:/GBM/results"

df = pd.read_csv(f"{ROOT}/receptor_tls_detection.csv")
rx = pd.read_excel("E:/GBM/ti tianran2.xlsx")

cat_names = ["Glutamate", "GABA/Gly", "Cholinergic", "DA/NE", "Serotonin"]
cat_colors = {"Glutamate": "#E64A19", "GABA/Gly": "#1B5E20",
              "Cholinergic": "#0D47A1", "DA/NE": "#6A1B9A", "Serotonin": "#BF360C"}

cat_map = {}
for i, col in enumerate(rx.columns):
    for g in rx[col].dropna():
        cat_map[g] = cat_names[i]

df["category"] = df["gene"].map(cat_map)
df = df.dropna(subset=["category"])

cols_data = []
max_n = 0
for cat in cat_names:
    sub = df[df["category"] == cat].sort_values("detect_rate", ascending=False)
    cols_data.append(sub)
    max_n = max(max_n, len(sub))

fig, ax = plt.subplots(figsize=(8, max_n * 0.20 + 0.6))

vmax = df["mean_expr"].quantile(0.95)
norm = Normalize(vmin=0, vmax=vmax)
cmap = mpl.colormaps["YlOrRd"]

for ci, (cat, sub) in enumerate(zip(cat_names, cols_data)):
    for i, (_, r) in enumerate(sub.iterrows()):
        s = r["detect_rate"] * 160 + 8  # size = detection rate
        c = cmap(norm(r["mean_expr"]))   # color = expression level
        ax.scatter(ci, i, s=s, c=[c], alpha=0.85, edgecolors="white", linewidths=0.3)
        ax.text(ci + 0.38, i, r["gene"], fontsize=5.5, va="center")

ax.set_xticks(range(5))
ax.set_xticklabels(cat_names, fontsize=7)
for i, c in enumerate(cat_names):
    ax.get_xticklabels()[i].set_color(cat_colors[c])
ax.tick_params(axis="x", length=0)
ax.set_yticks([])
ax.set_xlim(-0.8, 4.8)
ax.invert_yaxis()

# Colorbar for expression
cbar = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                     shrink=0.35, aspect=14, pad=0.02)
cbar.set_label("Mean expr.", fontsize=6)
cbar.ax.tick_params(labelsize=5.5)

# Size legend — placed left of first column, no overlap
for pct, label, xoff in [(25, "25%", -0.60), (50, "50%", -0.45), (75, "75%", -0.25)]:
    s = pct / 100 * 160 + 8
    ax.scatter(xoff, -1.2, s=s, c="grey", alpha=0.35,
               edgecolors="black", linewidths=0.3, clip_on=False)
    ax.text(xoff, -0.6, label, fontsize=5.5, ha="center")

fig.tight_layout()
for fmt in ["jpg", "svg", "pdf"]:
    fig.savefig(f"{ROOT}/fig_receptor_dotplot.{fmt}",
                dpi=600 if fmt == "jpg" else None, bbox_inches="tight")
plt.close()
print("Done: fig_receptor_dotplot.{jpg,svg,pdf}")
