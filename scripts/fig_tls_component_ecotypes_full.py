from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from sklearn.preprocessing import StandardScaler

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.6,
        "legend.frameon": False,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    }
)

ROOT = Path(r"E:/GBM/results")

ECO_ORDER = ["E1", "E2", "E3", "E4"]
ECO_LABELS = {
    "E1": "Lymphocyte TLS",
    "E2": "ILC-enriched TLS",
    "E3": "Myeloid-vascular TLS",
    "E4": "Glial-CD4 TLS",
}
ECO_SHORT = {
    "E1": "Lymphocyte",
    "E2": "ILC-enriched",
    "E3": "Myeloid-vascular",
    "E4": "Glial-CD4",
}
ECO_COLORS = {
    "E1": "#c44e52",
    "E2": "#4c72b0",
    "E3": "#55a868",
    "E4": "#dd8452",
}


def save_multi(fig, stem: str):
    fig.savefig(ROOT / f"{stem}.jpg", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def jitter(n, spread=0.14, seed=42):
    rng = np.random.RandomState(seed)
    return rng.uniform(-spread, spread, n)


def load_data():
    basis = pd.read_csv(ROOT / "tls_compnmf_rank4_basis.csv").set_index("cell_type")
    weights = pd.read_csv(ROOT / "tls_compnmf_rank4_unit_weights.csv")
    metrics = pd.read_csv(ROOT / "tls_compnmf_rank_metrics.csv")
    ann = pd.read_csv(ROOT / "tls_component_features_per_component.csv")
    meta = pd.read_csv(ROOT / "tls_pseudobulk_component_metadata.csv")
    summary = pd.read_csv(ROOT / "tls_compnmf_rank4_ecotype_summary.csv")

    ann = ann.merge(meta[["unit_id", "n_spots", "sample"]], on="unit_id", how="left")
    ann["sample"] = ann["sample_x"].fillna(ann["sample_y"]) if "sample_x" in ann.columns else ann["sample"]
    for col in ["sample_x", "sample_y"]:
        if col in ann.columns:
            ann = ann.drop(columns=col)
    ann["ecotype_name"] = ann["dominant_ecotype"].map(ECO_LABELS)
    ann["dominant_ecotype"] = pd.Categorical(ann["dominant_ecotype"], categories=ECO_ORDER, ordered=True)
    weights["dominant_ecotype"] = pd.Categorical(weights["dominant_ecotype"], categories=ECO_ORDER, ordered=True)
    return basis, weights, metrics, ann, summary


def plot_nmf_basis(basis, summary):
    basis = basis.loc[
        [
            "B",
            "Plasma",
            "CD4_T",
            "CD8_T",
            "NK",
            "Dendritic",
            "ILC1",
            "ILC2",
            "ILC3",
            "Macrophage",
            "Vascular",
            "Glial",
            "Glioma",
        ]
    ]
    basis_z = StandardScaler().fit_transform(basis.values.T).T
    vmax = np.abs(basis_z).max()

    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(basis_z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(basis.index)))
    ax.set_yticklabels(basis.index, fontsize=8)
    ax.set_xticks(range(len(ECO_ORDER)))
    ax.set_xticklabels([ECO_SHORT[e] for e in ECO_ORDER], rotation=25, ha="right", fontsize=8)
    for idx, eco in enumerate(ECO_ORDER):
        ax.get_xticklabels()[idx].set_color(ECO_COLORS[eco])
        n_units = int(summary.loc[summary["ecotype"] == eco, "n_units"].iloc[0])
        ax.text(idx, len(basis.index) + 0.45, f"n={n_units}", ha="center", va="top", fontsize=7, color="#555555")
    ax.set_title("NMF Basis Heatmap", fontsize=10, pad=8)
    cbar = plt.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
    cbar.set_label("Row z-score", fontsize=8)
    ax.set_ylim(len(basis.index) + 0.9, -0.5)
    fig.tight_layout()
    save_multi(fig, "fig_nmf_basis")


def plot_nmf_coefficients(weights):
    eco_cols = ECO_ORDER
    df = weights.copy()
    df["dominant_rank"] = df["dominant_ecotype"].cat.codes
    df = df.sort_values(["dominant_rank", "dominant_weight"], ascending=[True, False]).reset_index(drop=True)
    H = df[eco_cols].to_numpy().T

    counts = [(df["dominant_ecotype"] == eco).sum() for eco in ECO_ORDER]

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    im = ax.imshow(H, aspect="auto", cmap="YlOrRd", vmin=0, vmax=np.quantile(H, 0.98))
    ax.set_yticks(range(len(ECO_ORDER)))
    ax.set_yticklabels([f"{ECO_SHORT[e]} (n={n})" for e, n in zip(ECO_ORDER, counts)], fontsize=8)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(df)} TLS components", fontsize=9)
    ax.set_title("NMF Coefficients by Component", fontsize=10, pad=8)

    cursor = 0
    for eco, n in zip(ECO_ORDER, counts):
        ax.axvline(cursor - 0.5, color="white", linewidth=0.8)
        ax.add_patch(
            plt.Rectangle((cursor - 0.5, -0.88), n, 0.22, color=ECO_COLORS[eco], clip_on=False)
        )
        ax.text(cursor + n / 2 - 0.5, -1.1, ECO_SHORT[eco], ha="center", va="top", fontsize=7, color=ECO_COLORS[eco])
        cursor += n
    ax.axvline(cursor - 0.5, color="white", linewidth=0.8)

    cbar = plt.colorbar(im, ax=ax, shrink=0.65, pad=0.02)
    cbar.set_label("Normalized NMF weight", fontsize=8)
    fig.tight_layout()
    save_multi(fig, "fig_nmf_coefficients")


def plot_rank_selection(metrics):
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.plot(metrics["rank"], metrics["cophenetic"], marker="o", color="#1f77b4", label="Cophenetic")
    ax.plot(metrics["rank"], metrics["silhouette_consensus"], marker="s", color="#d62728", label="Consensus silhouette")
    ax.plot(metrics["rank"], metrics["dispersion"], marker="^", color="#2ca02c", label="Dispersion")
    ax.axvline(4, linestyle="--", linewidth=0.8, color="#666666")
    ax.text(4.04, metrics[["cophenetic", "silhouette_consensus", "dispersion"]].max().max(), "K=4", fontsize=8, va="top")
    ax.set_xticks(metrics["rank"])
    ax.set_xlabel("NMF rank (K)", fontsize=9)
    ax.set_ylabel("Metric value", fontsize=9)
    ax.set_title("Rank Selection", fontsize=10, pad=8)
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    save_multi(fig, "fig_nmf_rank_selection")


def box_scatter(ax, data_groups, colors, labels, ylabel):
    bp = ax.boxplot(
        data_groups,
        patch_artist=True,
        widths=0.58,
        medianprops={"color": "black", "linewidth": 0.8},
        flierprops={"marker": "none"},
        whiskerprops={"linewidth": 0.6},
        capprops={"linewidth": 0.6},
    )
    for idx, (patch, vals) in enumerate(zip(bp["boxes"], data_groups)):
        patch.set_facecolor(colors[idx])
        patch.set_alpha(0.35)
        x = np.full(len(vals), idx + 1) + jitter(len(vals), seed=42 + idx)
        ax.scatter(x, vals, s=6, alpha=0.24, color=colors[idx], linewidth=0)
    ax.set_xticklabels(labels, rotation=24, ha="right", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)


def plot_ecotype_features(ann):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8))
    labels = [ECO_SHORT[e] for e in ECO_ORDER]
    colors = [ECO_COLORS[e] for e in ECO_ORDER]

    maturity = [ann.loc[ann["dominant_ecotype"] == e, "maturity_score"].dropna().values for e in ECO_ORDER]
    ilc_total = [ann.loc[ann["dominant_ecotype"] == e, "ILC_total"].dropna().values for e in ECO_ORDER]
    n_spots = [ann.loc[ann["dominant_ecotype"] == e, "n_spots"].dropna().values for e in ECO_ORDER]

    box_scatter(axes[0], maturity, colors, labels, "Maturity score")
    box_scatter(axes[1], ilc_total, colors, labels, "ILC total abundance")
    box_scatter(axes[2], n_spots, colors, labels, "Spots per component")
    axes[0].set_title("Ecotype Features", fontsize=10, pad=8)
    fig.tight_layout()
    save_multi(fig, "fig_ecotype_features")


def plot_maturity_vs_ilc(ann):
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for eco in ECO_ORDER:
        sub = ann[ann["dominant_ecotype"] == eco]
        ax.scatter(
            sub["maturity_score"],
            sub["ILC_total"],
            s=10,
            alpha=0.35,
            color=ECO_COLORS[eco],
            label=f"{ECO_SHORT[eco]} (n={len(sub)})",
            linewidth=0,
        )
    ax.set_xlabel("Maturity score", fontsize=9)
    ax.set_ylabel("ILC total abundance", fontsize=9)
    ax.set_title("Maturity vs ILC", fontsize=10, pad=8)
    ax.legend(fontsize=7, loc="lower right", markerscale=1.3)
    fig.tight_layout()
    save_multi(fig, "fig_maturity_vs_ilc")


def plot_ilc_proportions(ann):
    frac = ann.copy()
    frac["ILC1_frac"] = frac["ILC1"] / frac["ILC_total"].clip(lower=1e-8)
    frac["ILC2_frac"] = frac["ILC2"] / frac["ILC_total"].clip(lower=1e-8)
    frac["ILC3_frac"] = frac["ILC3"] / frac["ILC_total"].clip(lower=1e-8)

    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.6), sharey=True)
    for ax, feat, title, fill in zip(
        axes,
        ["ILC1_frac", "ILC2_frac", "ILC3_frac"],
        ["ILC1 fraction", "ILC2 fraction", "ILC3 fraction"],
        ["#7a9e9f", "#c17c74", "#7b8fce"],
    ):
        data = [frac.loc[frac["dominant_ecotype"] == eco, feat].dropna().values for eco in ECO_ORDER]
        bp = ax.boxplot(
            data,
            patch_artist=True,
            widths=0.58,
            medianprops={"color": "black", "linewidth": 0.8},
            flierprops={"marker": "none"},
        )
        for idx, (patch, vals) in enumerate(zip(bp["boxes"], data)):
            patch.set_facecolor(fill)
            patch.set_alpha(0.4)
            x = np.full(len(vals), idx + 1) + jitter(len(vals), seed=80 + idx)
            ax.scatter(x, vals, s=5, alpha=0.2, color=ECO_COLORS[ECO_ORDER[idx]], linewidth=0)
        ax.set_xticklabels([ECO_SHORT[e] for e in ECO_ORDER], rotation=24, ha="right", fontsize=7)
        ax.set_title(title, fontsize=9, pad=6)
    axes[0].set_ylabel("Fraction within ILC total", fontsize=8)
    fig.tight_layout()
    save_multi(fig, "fig_ilc_proportions")


def plot_ecotype_counts(ann):
    comp_counts = ann.groupby("dominant_ecotype").size().reindex(ECO_ORDER)
    sample_counts = ann.groupby("dominant_ecotype")["sample"].nunique().reindex(ECO_ORDER)
    x = np.arange(len(ECO_ORDER))
    width = 0.34

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.bar(x - width / 2, comp_counts.values, width=width, color=[ECO_COLORS[e] for e in ECO_ORDER], alpha=0.9, label="Components")
    ax.bar(x + width / 2, sample_counts.values, width=width, color="#999999", alpha=0.8, label="Samples")
    ax.set_xticks(x)
    ax.set_xticklabels([ECO_SHORT[e] for e in ECO_ORDER], rotation=22, ha="right", fontsize=8)
    ax.set_ylabel("Count", fontsize=9)
    ax.set_title("Ecotype Coverage", fontsize=10, pad=8)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    save_multi(fig, "fig_ecotype_counts")


def plot_receptor_dotplot():
    rx = pd.read_excel(r"E:/GBM/ti tianran2.xlsx")
    rx_map = {}
    cat_names = ["Glutamate", "GABA/Gly", "Cholinergic", "DA/NE", "Serotonin"]
    cat_colors = ["#c95f37", "#4f7f3d", "#3e6aa7", "#8b5a9a", "#b55a2a"]
    for cat, col in zip(cat_names, rx.columns):
        for gene in rx[col].dropna():
            rx_map[gene] = cat

    df = pd.read_csv(ROOT / "receptor_pseudobulk_ecotype.csv")
    df["category"] = df["gene"].map(rx_map)
    df = df.dropna(subset=["category"])
    df["ecotype"] = pd.Categorical(df["ecotype"], categories=ECO_ORDER, ordered=True)

    groups = []
    max_n = 0
    for cat in cat_names:
        for eco in ECO_ORDER:
            sub = df[(df["category"] == cat) & (df["ecotype"] == eco)].sort_values("detect_rate", ascending=False)
            groups.append((cat, eco, sub))
            max_n = max(max_n, len(sub))

    fig, ax = plt.subplots(figsize=(15.5, max_n * 0.19 + 1.6))
    vmax = df["mean_expr"].quantile(0.95)
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = mpl.colormaps["YlOrRd"]

    for ci, (cat, eco, sub) in enumerate(groups):
        first_of_cat = ci % 4 == 0
        for ri, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(
                ci,
                ri,
                s=row["detect_rate"] * 145 + 8,
                c=[cmap(norm(row["mean_expr"]))],
                alpha=0.88,
                edgecolors="white",
                linewidths=0.3,
            )
            if first_of_cat:
                ax.text(ci + 0.38, ri, row["gene"], fontsize=5, va="center")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([ECO_SHORT[eco] for _, eco, _ in groups], fontsize=5.6, rotation=28, ha="right")
    for idx, cat in enumerate(cat_names):
        mid = idx * 4 + 1.5
        ax.text(mid, -1.45, cat, fontsize=7, fontweight="bold", color=cat_colors[idx], ha="center", va="bottom")
        if idx > 0:
            ax.axvline(idx * 4 - 0.5, color="#d9d9d9", linewidth=0.7, linestyle="--")
    ax.set_yticks([])
    ax.set_xlim(-0.8, len(groups) - 0.2)
    ax.invert_yaxis()
    ax.set_title("Receptor Programs Across TLS Ecotypes", fontsize=10, pad=18)

    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.28, aspect=18, pad=0.02)
    cbar.set_label("Mean CPM", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    for pct, xoff in zip([25, 50, 75], [-0.65, -0.48, -0.28]):
        ax.scatter(xoff, -0.95, s=pct / 100 * 145 + 8, c="gray", alpha=0.35, edgecolors="black", linewidths=0.3, clip_on=False)
        ax.text(xoff, -0.18, f"{pct}%", fontsize=5.5, ha="center")

    fig.tight_layout()
    save_multi(fig, "fig_receptor_ecotype_dotplot")


def main():
    basis, weights, metrics, ann, summary = load_data()
    plot_nmf_basis(basis, summary)
    plot_nmf_coefficients(weights)
    plot_rank_selection(metrics)
    plot_ecotype_features(ann)
    plot_maturity_vs_ilc(ann)
    plot_ilc_proportions(ann)
    plot_ecotype_counts(ann)
    plot_receptor_dotplot()
    print("Saved TLS component ecotype figure set.")


if __name__ == "__main__":
    main()
