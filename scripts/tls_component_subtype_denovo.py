"""De novo TLS component subtype analysis for GBM.

Rebuilds TLS component states from scratch using:
1. SpaLinker-defined TLS spots
2. component extraction by spatial connectivity
3. component composition clustering on CLR-transformed cell2location composition
4. maturity and ILC enrichment summaries
"""

from __future__ import annotations

import warnings
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import RobustScaler

warnings.filterwarnings("ignore")

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.5,
        "legend.frameon": False,
    }
)

H5AD = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata_consolidated")
TLS_DIR = Path(r"E:/GBM/results/tls_consolidated")
OUT = Path(r"E:/GBM/results")

SKIP = {"GSE194329_DMG1", "GSE194329_DMG2", "GSE194329_DMG3", "GSE194329_DMG4", "GSE194329_DMG5"}
MIN_SHARED = 100
MIN_TLS = 5
MIN_COMPONENT = 5
K_NEIGH = 7
K_RANGE = range(2, 9)


def to_numpy(x) -> np.ndarray:
    if hasattr(x, "values"):
        x = x.values
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x)


def clr_transform(prop: np.ndarray) -> np.ndarray:
    safe = np.clip(prop, 1e-8, None)
    logp = np.log(safe)
    return logp - logp.mean(axis=1, keepdims=True)


def collect_components() -> tuple[pd.DataFrame, list[str]]:
    rows = []
    ct_names: list[str] | None = None

    for h5_path in sorted(H5AD.glob("*.h5ad")):
        if h5_path.stem in SKIP:
            continue

        tls_csv = TLS_DIR / h5_path.stem / "tls_spot_scores_official_relaxed.csv"
        if not tls_csv.exists():
            continue

        tls = pd.read_csv(tls_csv)
        if "barcode" in tls.columns:
            tls = tls.set_index("barcode")

        adata = ad.read_h5ad(h5_path)
        shared = adata.obs_names.intersection(tls.index)
        if len(shared) < MIN_SHARED:
            continue

        adata = adata[shared].copy()
        tls = tls.loc[shared]
        tls_mask = (tls["TLS.region"] == "TLS").to_numpy()
        if int(tls_mask.sum()) < MIN_TLS:
            continue

        q05 = to_numpy(adata.obsm["c2l_ilc_q05"]).astype(np.float64)
        prop = q05 / (q05.sum(axis=1, keepdims=True) + 1e-8)
        coords = to_numpy(adata.obsm["spatial"]).astype(np.float64)
        if ct_names is None:
            ct_names = list(adata.uns["c2l_ilc_cell_types"])

        tls_idx = np.where(tls_mask)[0]
        tls_coords = coords[tls_idx]
        k = min(K_NEIGH, len(tls_idx))
        adj = kneighbors_graph(tls_coords, n_neighbors=k, mode="connectivity", include_self=True)
        n_comp, labels = connected_components(adj, directed=False)

        for comp_id in range(n_comp):
            comp_local = labels == comp_id
            comp_idx = tls_idx[comp_local]
            if len(comp_idx) < MIN_COMPONENT:
                continue

            comp_q05 = q05[comp_idx]
            comp_prop = prop[comp_idx]
            row = {
                "sample": h5_path.stem,
                "component_id": comp_id,
                "tls_size": len(comp_idx),
                "tls_score_mean": float(tls.iloc[comp_idx]["TLS.score"].mean()) if "TLS.score" in tls.columns else np.nan,
            }
            for i, ct in enumerate(ct_names):
                row[f"{ct}_mean"] = float(comp_q05[:, i].mean())
                row[f"{ct}_prop"] = float(comp_prop[:, i].mean())
            rows.append(row)

    if ct_names is None:
        raise RuntimeError("No eligible TLS components found.")

    return pd.DataFrame(rows), ct_names


def cluster_components(df: pd.DataFrame, ct_names: list[str]) -> tuple[pd.DataFrame, int, pd.DataFrame]:
    prop_cols = [f"{ct}_prop" for ct in ct_names]
    X = df[prop_cols].to_numpy()
    X_clr = clr_transform(X)
    X_scaled = RobustScaler().fit_transform(X_clr)

    scores = []
    best_k = None
    best_score = -np.inf
    best_labels = None
    for k in K_RANGE:
        labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X_scaled)
        score = silhouette_score(X_scaled, labels)
        scores.append({"k": k, "silhouette": score})
        if score > best_score:
            best_k = k
            best_score = score
            best_labels = labels

    out = df.copy()
    out["subtype"] = best_labels

    score_df = pd.DataFrame(scores)
    return out, int(best_k), score_df


def summarize_components(df: pd.DataFrame, ct_names: list[str]) -> pd.DataFrame:
    comp = df.copy()
    comp["ILC_total"] = comp["ILC1_mean"] + comp["ILC2_mean"] + comp["ILC3_mean"]
    comp["maturity_score"] = (
        comp["B_mean"]
        + comp["Plasma_mean"]
        + comp["CD4_T_mean"]
        + comp["CD8_T_mean"]
        + comp["Dendritic_mean"]
        + comp["NK_mean"]
    )
    comp["ILC1_frac"] = comp["ILC1_mean"] / (comp["ILC_total"] + 1e-8)
    comp["ILC2_frac"] = comp["ILC2_mean"] / (comp["ILC_total"] + 1e-8)
    comp["ILC3_frac"] = comp["ILC3_mean"] / (comp["ILC_total"] + 1e-8)

    rows = []
    for subtype, sub in comp.groupby("subtype"):
        mean_delta = []
        for ct in ct_names:
            mean_delta.append(sub[f"{ct}_prop"].mean())
        order = np.argsort(mean_delta)[::-1]
        rows.append(
            {
                "subtype": subtype,
                "n_components": len(sub),
                "n_samples": sub["sample"].nunique(),
                "tls_size_median": sub["tls_size"].median(),
                "tls_score_mean": sub["tls_score_mean"].mean(),
                "maturity_score_mean": sub["maturity_score"].mean(),
                "ILC_total_mean": sub["ILC_total"].mean(),
                "ILC1_frac_median": sub["ILC1_frac"].median(),
                "ILC2_frac_median": sub["ILC2_frac"].median(),
                "ILC3_frac_median": sub["ILC3_frac"].median(),
                "top_1": ct_names[order[0]],
                "top_2": ct_names[order[1]],
                "top_3": ct_names[order[2]],
                **{ct: sub[f"{ct}_prop"].mean() for ct in ct_names},
            }
        )

    return comp, pd.DataFrame(rows).sort_values("maturity_score_mean", ascending=False).reset_index(drop=True)


def plot_outputs(summary: pd.DataFrame, ct_names: list[str], score_df: pd.DataFrame, best_k: int) -> None:
    fig = plt.figure(figsize=(13, 8))

    ax1 = fig.add_axes([0.05, 0.52, 0.52, 0.40])
    mat = summary[ct_names].to_numpy().T
    vmax = max(abs(mat.min()), abs(mat.max()), 0.05)
    im = ax1.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0, vmax=mat.max())
    ax1.set_yticks(range(len(ct_names)))
    ax1.set_yticklabels(ct_names, fontsize=6)
    ax1.set_xticks(range(len(summary)))
    ax1.set_xticklabels(
        [f"S{int(summary.loc[i, 'subtype'])}\n{summary.loc[i, 'top_1']}\nN={int(summary.loc[i, 'n_components'])}" for i in range(len(summary))],
        fontsize=6,
    )
    ax1.set_title(
        f"a  De novo TLS component subtype composition (best K={best_k})",
        fontsize=9,
        fontweight="bold",
        loc="left",
    )
    plt.colorbar(im, ax=ax1, shrink=0.7).set_label("mean component proportion", fontsize=6)

    ax2 = fig.add_axes([0.64, 0.52, 0.13, 0.40])
    ax2.bar(range(len(summary)), summary["maturity_score_mean"], color="#c44e52", edgecolor="black", linewidth=0.3)
    ax2.set_xticks(range(len(summary)))
    ax2.set_xticklabels([f"S{int(x)}" for x in summary["subtype"]], rotation=90, fontsize=6)
    ax2.set_ylabel("maturity score", fontsize=7)
    ax2.set_title("b  Maturity", fontsize=9, fontweight="bold", loc="left")

    ax3 = fig.add_axes([0.82, 0.52, 0.13, 0.40])
    frac_mat = summary[["ILC1_frac_median", "ILC2_frac_median", "ILC3_frac_median"]].to_numpy().T
    im3 = ax3.imshow(frac_mat, aspect="auto", cmap="YlGnBu", vmin=0.28, vmax=0.38)
    ax3.set_yticks(range(3))
    ax3.set_yticklabels(["ILC1", "ILC2", "ILC3"], fontsize=6)
    ax3.set_xticks(range(len(summary)))
    ax3.set_xticklabels([f"S{int(x)}" for x in summary["subtype"]], rotation=90, fontsize=6)
    ax3.set_title("c  ILC subtype fraction", fontsize=9, fontweight="bold", loc="left")
    plt.colorbar(im3, ax=ax3, shrink=0.7)

    ax4 = fig.add_axes([0.08, 0.10, 0.24, 0.25])
    ax4.plot(score_df["k"], score_df["silhouette"], marker="o", color="black")
    ax4.axvline(best_k, color="#c44e52", linestyle="--", linewidth=1)
    ax4.set_xlabel("K", fontsize=7)
    ax4.set_ylabel("silhouette", fontsize=7)
    ax4.set_title("d  Model selection", fontsize=9, fontweight="bold", loc="left")

    ax5 = fig.add_axes([0.40, 0.10, 0.24, 0.25])
    ax5.scatter(summary["maturity_score_mean"], summary["ILC_total_mean"], s=36, color="#4c72b0")
    for _, r in summary.iterrows():
        ax5.text(r["maturity_score_mean"], r["ILC_total_mean"], f"S{int(r['subtype'])}", fontsize=7)
    ax5.set_xlabel("maturity score", fontsize=7)
    ax5.set_ylabel("ILC total", fontsize=7)
    ax5.set_title("e  Maturity vs ILC", fontsize=9, fontweight="bold", loc="left")

    ax6 = fig.add_axes([0.72, 0.10, 0.23, 0.25])
    ax6.bar(range(len(summary)), summary["n_samples"], color="#55a868", edgecolor="black", linewidth=0.3)
    ax6.set_xticks(range(len(summary)))
    ax6.set_xticklabels([f"S{int(x)}" for x in summary["subtype"]], rotation=90, fontsize=6)
    ax6.set_ylabel("samples", fontsize=7)
    ax6.set_title("f  Sample support", fontsize=9, fontweight="bold", loc="left")

    fig.savefig(OUT / "fig_tls_component_subtype_denovo.jpg", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "fig_tls_component_subtype_denovo.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_tls_component_subtype_denovo.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df, ct_names = collect_components()
    clustered, best_k, score_df = cluster_components(df, ct_names)
    comp, summary = summarize_components(clustered, ct_names)

    comp.to_csv(OUT / "tls_components_denovo.csv", index=False)
    summary.to_csv(OUT / "tls_component_subtype_denovo_summary.csv", index=False)
    score_df.to_csv(OUT / "tls_component_subtype_denovo_model_selection.csv", index=False)
    plot_outputs(summary, ct_names, score_df, best_k)

    print(f"components={len(comp)}")
    print(f"best_k={best_k}")
    print(summary[["subtype", "n_components", "n_samples", "maturity_score_mean", "ILC_total_mean", "ILC1_frac_median", "ILC2_frac_median", "ILC3_frac_median", "top_1", "top_2", "top_3"]].to_string(index=False))


if __name__ == "__main__":
    main()
