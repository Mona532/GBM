from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree


def feature_mask(names: list[str], keywords: list[str]) -> np.ndarray:
    lowered = [name.lower() for name in names]
    return np.array([any(key.lower() in name for key in keywords) for name in lowered], dtype=bool)


def get_dense_columns(matrix: sparse.spmatrix | np.ndarray, mask: np.ndarray) -> np.ndarray:
    if mask.sum() == 0:
        n_obs = matrix.shape[0]
        return np.zeros((n_obs, 0), dtype=float)
    if sparse.issparse(matrix):
        return matrix[:, mask].toarray()
    return np.asarray(matrix[:, mask])


def robust_z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        sd = np.nanstd(x)
        if sd == 0 or np.isnan(sd):
            return np.zeros_like(x)
        return (x - np.nanmean(x)) / sd
    return (x - med) / (1.4826 * mad)


def smooth_knn(values: np.ndarray, coords: np.ndarray, k: int = 6) -> np.ndarray:
    n = len(values)
    if n == 0:
        return values
    k_eff = min(k + 1, n)
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k_eff)
    if k_eff == 1:
        idx = idx[:, None]
    return values[idx].mean(axis=1)


def connected_hotspots(coords: np.ndarray, hot: np.ndarray, k: int = 6) -> np.ndarray:
    n = len(hot)
    labels = np.zeros(n, dtype=int)
    hot_idx = np.flatnonzero(hot)
    if hot_idx.size == 0:
        return labels
    tree = cKDTree(coords)
    _, nn = tree.query(coords, k=min(k + 1, n))
    if nn.ndim == 1:
        nn = nn[:, None]
    current = 0
    hot_set = set(hot_idx.tolist())
    unseen = set(hot_idx.tolist())
    while unseen:
        seed = unseen.pop()
        current += 1
        stack = [seed]
        labels[seed] = current
        while stack:
            node = stack.pop()
            for nbr in nn[node]:
                nbr = int(nbr)
                if nbr in hot_set and labels[nbr] == 0:
                    labels[nbr] = current
                    if nbr in unseen:
                        unseen.remove(nbr)
                    stack.append(nbr)
    return labels


def plot_tls(coords: np.ndarray, tls_score: np.ndarray, hot: np.ndarray, output_png: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
    sc = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=tls_score,
        s=10,
        cmap="inferno",
        linewidths=0,
    )
    if hot.any():
        ax.scatter(
            coords[hot, 0],
            coords[hot, 1],
            facecolors="none",
            edgecolors="#39a0ed",
            s=20,
            linewidths=0.6,
        )
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("TLS_score_relaxed")
    fig.tight_layout()
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)


def score_sample(path: Path, out_dir: Path, top_quantile: float, hotspot_quantile: float, k: int) -> dict:
    adata = ad.read_h5ad(path)
    var_names = adata.var_names.astype(str).tolist()
    feature_types = adata.var["feature_types"].astype(str).tolist()
    x = adata.X

    cell_state_mask = np.array([ft == "Cell state abundances" for ft in feature_types], dtype=bool)
    niche_mask = np.array([ft == "Spatial niche abundances" for ft in feature_types], dtype=bool)
    histo_mask = np.array([ft == "Histopath annotation overlap" for ft in feature_types], dtype=bool)

    cell_state_names = [name for name, keep in zip(var_names, cell_state_mask) if keep]
    niche_names = [name for name, keep in zip(var_names, niche_mask) if keep]
    histo_names = [name for name, keep in zip(var_names, histo_mask) if keep]

    b_mask_local = feature_mask(cell_state_names, ["b cells", "plasma"])
    t_mask_local = feature_mask(cell_state_names, ["cd4", "cd8", "t cell"])
    dc_mask_local = feature_mask(cell_state_names, ["dendritic"])
    niche_mask_local = feature_mask(niche_names, ["immune"])
    vascular_mask_local = feature_mask(histo_names, ["hyperplastic blood vessels", "microvascular proliferation"])
    necrosis_mask_local = feature_mask(histo_names, ["necrosis", "perinecrotic", "pseudopalisading"])

    cell_state_values = get_dense_columns(x, cell_state_mask)
    niche_values = get_dense_columns(x, niche_mask)
    histo_values = get_dense_columns(x, histo_mask)

    b_plasma = cell_state_values[:, b_mask_local].sum(axis=1) if b_mask_local.any() else np.zeros(adata.n_obs)
    t_cell = cell_state_values[:, t_mask_local].sum(axis=1) if t_mask_local.any() else np.zeros(adata.n_obs)
    dendritic = cell_state_values[:, dc_mask_local].sum(axis=1) if dc_mask_local.any() else np.zeros(adata.n_obs)
    immune_niche = niche_values[:, niche_mask_local].sum(axis=1) if niche_mask_local.any() else np.zeros(adata.n_obs)
    vascular = histo_values[:, vascular_mask_local].sum(axis=1) if vascular_mask_local.any() else np.zeros(adata.n_obs)
    necrosis = histo_values[:, necrosis_mask_local].sum(axis=1) if necrosis_mask_local.any() else np.zeros(adata.n_obs)

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    b_smooth = smooth_knn(b_plasma, coords, k=k)
    t_smooth = smooth_knn(t_cell, coords, k=k)
    dc_smooth = smooth_knn(dendritic, coords, k=k)
    immune_smooth = smooth_knn(immune_niche, coords, k=k)

    bt_coloc = np.sqrt(np.maximum(b_smooth, 0) * np.maximum(t_smooth, 0))
    btd_coloc = np.cbrt(np.maximum(b_smooth, 0) * np.maximum(t_smooth, 0) * np.maximum(dc_smooth + 1e-8, 0))

    score_components = pd.DataFrame(
        {
            "b_plasma_z": robust_z(np.log1p(b_smooth)),
            "t_cell_z": robust_z(np.log1p(t_smooth)),
            "dendritic_z": robust_z(np.log1p(dc_smooth)),
            "bt_coloc_z": robust_z(np.log1p(bt_coloc)),
            "btd_coloc_z": robust_z(np.log1p(btd_coloc)),
            "immune_niche_z": robust_z(np.log1p(immune_smooth)),
            "vascular_z": robust_z(np.log1p(vascular)) if vascular.sum() > 0 else np.zeros(adata.n_obs),
            "necrosis_z": robust_z(np.log1p(necrosis)) if necrosis.sum() > 0 else np.zeros(adata.n_obs),
        }
    )

    tls_score = (
        0.25 * score_components["b_plasma_z"].to_numpy()
        + 0.25 * score_components["t_cell_z"].to_numpy()
        + 0.10 * score_components["dendritic_z"].to_numpy()
        + 0.20 * score_components["bt_coloc_z"].to_numpy()
        + 0.10 * score_components["btd_coloc_z"].to_numpy()
        + 0.10 * score_components["immune_niche_z"].to_numpy()
        + 0.05 * score_components["vascular_z"].to_numpy()
        - 0.05 * score_components["necrosis_z"].to_numpy()
    )

    eligible = (b_plasma > 0) | np.array([False] * adata.n_obs)
    positive_t = t_cell > 0
    score_cut = np.quantile(tls_score[eligible & positive_t], top_quantile) if np.any(eligible & positive_t) else np.inf
    coloc_cut = np.quantile(bt_coloc[eligible & positive_t], hotspot_quantile) if np.any(eligible & positive_t) else np.inf
    hotspot = (tls_score >= score_cut) & (bt_coloc >= coloc_cut) & eligible & positive_t
    hotspot_component = connected_hotspots(coords, hotspot, k=k)

    sample_dir = out_dir / path.stem
    sample_dir.mkdir(parents=True, exist_ok=True)

    spot_df = pd.DataFrame(
        {
            "barcode": adata.obs_names.astype(str),
            "sample_id": path.stem,
            "array_row": adata.obs["array_row"].to_numpy() if "array_row" in adata.obs else np.nan,
            "array_col": adata.obs["array_col"].to_numpy() if "array_col" in adata.obs else np.nan,
            "in_tissue": adata.obs["in_tissue"].to_numpy() if "in_tissue" in adata.obs else np.nan,
            "x": coords[:, 0],
            "y": coords[:, 1],
            "B_plasma_raw": b_plasma,
            "T_cell_raw": t_cell,
            "Dendritic_raw": dendritic,
            "Immune_niche_raw": immune_niche,
            "Vascular_histopath_raw": vascular,
            "Necrosis_histopath_raw": necrosis,
            "B_plasma_smooth": b_smooth,
            "T_cell_smooth": t_smooth,
            "Dendritic_smooth": dc_smooth,
            "BT_coloc": bt_coloc,
            "BTD_coloc": btd_coloc,
            "TLS_score_relaxed": tls_score,
            "TLS_hotspot_relaxed": hotspot,
            "TLS_hotspot_component": hotspot_component,
        }
    )
    spot_df = pd.concat([spot_df, score_components], axis=1)
    spot_df.to_csv(sample_dir / "tls_spot_scores.csv", index=False)

    manifest = {
        "sample_id": path.stem,
        "source_h5ad": str(path),
        "n_spots": int(adata.n_obs),
        "selected_b_plasma_features": [name for name, keep in zip(cell_state_names, b_mask_local) if keep],
        "selected_t_cell_features": [name for name, keep in zip(cell_state_names, t_mask_local) if keep],
        "selected_dendritic_features": [name for name, keep in zip(cell_state_names, dc_mask_local) if keep],
        "selected_immune_niche_features": [name for name, keep in zip(niche_names, niche_mask_local) if keep],
        "selected_vascular_histopath_features": [name for name, keep in zip(histo_names, vascular_mask_local) if keep],
        "selected_necrosis_histopath_features": [name for name, keep in zip(histo_names, necrosis_mask_local) if keep],
        "score_quantile_cutoff": float(score_cut) if np.isfinite(score_cut) else None,
        "coloc_quantile_cutoff": float(coloc_cut) if np.isfinite(coloc_cut) else None,
    }
    with open(sample_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)

    plot_tls(coords, tls_score, hotspot, sample_dir / "tls_score_map.png", path.stem)

    return {
        "sample_id": path.stem,
        "n_spots": int(adata.n_obs),
        "has_b_or_plasma": bool((b_plasma > 0).any()),
        "has_t_cell": bool((t_cell > 0).any()),
        "n_hotspots": int(hotspot.sum()),
        "n_hotspot_components": int(hotspot_component.max()),
        "tls_score_mean": float(np.mean(tls_score)),
        "tls_score_median": float(np.median(tls_score)),
        "tls_score_q95": float(np.quantile(tls_score, 0.95)),
        "bt_coloc_q95": float(np.quantile(bt_coloc, 0.95)),
        "b_plasma_features": "; ".join([name for name, keep in zip(cell_state_names, b_mask_local) if keep]),
        "t_cell_features": "; ".join([name for name, keep in zip(cell_state_names, t_mask_local) if keep]),
        "dendritic_features": "; ".join([name for name, keep in zip(cell_state_names, dc_mask_local) if keep]),
        "immune_niche_features": "; ".join([name for name, keep in zip(niche_names, niche_mask_local) if keep]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Relaxed TLS scoring for Visium h5ad files.")
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--top-quantile", type=float, default=0.95)
    parser.add_argument("--hotspot-quantile", type=float, default=0.75)
    parser.add_argument("--k", type=int, default=6)
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("*.h5ad"))
    if not files:
        raise SystemExit(f"No .h5ad files found under {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for path in files:
        print(f"[score] {path.name}")
        summaries.append(
            score_sample(
                path=path,
                out_dir=args.output_dir,
                top_quantile=args.top_quantile,
                hotspot_quantile=args.hotspot_quantile,
                k=args.k,
            )
        )

    summary_df = pd.DataFrame(summaries).sort_values("sample_id")
    summary_df.to_csv(args.output_dir / "tls_relaxed_summary.csv", index=False)
    print(f"[done] {summary_df.shape[0]} samples written to {args.output_dir}")


if __name__ == "__main__":
    main()
