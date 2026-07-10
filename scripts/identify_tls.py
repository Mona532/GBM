from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph


# Cell types defining TLS immune signature
TLS_CELL_TYPES = ["B", "CD8_T", "Dendritic", "NK", "ILC1", "ILC2", "ILC3"]

# Neighbors per spot (Visium hexagonal grid, 6 neighbors + self)
K_NEIGHBORS = 7


def compute_tls_score(q05: np.ndarray, cell_types: np.ndarray) -> pd.DataFrame:
    """Compute per-spot TLS scores from q05 abundance matrix."""
    tls_idx = [list(cell_types).index(c) for c in TLS_CELL_TYPES]
    tls_abundance = q05[:, tls_idx]  # (n_spots, n_tls_types)

    spot_df = pd.DataFrame(
        tls_abundance,
        columns=[f"tls_{c}" for c in TLS_CELL_TYPES],
    )
    spot_df["tls_sum"] = tls_abundance.sum(axis=1)
    # Z-score within sample
    spot_df["tls_zscore"] = (spot_df["tls_sum"] - spot_df["tls_sum"].mean()) / spot_df["tls_sum"].std()
    # Proportion of total per spot
    total_per_spot = q05.sum(axis=1)
    spot_df["tls_ratio"] = spot_df["tls_sum"] / np.where(total_per_spot == 0, 1, total_per_spot)
    return spot_df


def build_spatial_graph(positions: np.ndarray, k: int = K_NEIGHBORS) -> csr_matrix:
    """Build k-nearest-neighbor spatial graph from 2D coordinates."""
    return kneighbors_graph(positions, n_neighbors=k, mode="connectivity", include_self=True)


def label_tls_spots(
    spot_df: pd.DataFrame,
    adjacency: csr_matrix,
    z_threshold: float = 1.0,
    min_cluster_size: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Identify TLS-positive spots and cluster them into individual TLS structures.

    Returns (tls_label, tls_cluster_id) arrays.
        tls_label: 1 = TLS-positive, 0 = not TLS
        tls_cluster_id: 0 = no TLS, >=1 = TLS cluster ID
    """
    n_spots = len(spot_df)
    # Step 1: candidate spots with high immune z-score
    candidates = spot_df["tls_zscore"].values >= z_threshold
    print(f"  candidates (z>={z_threshold}): {candidates.sum()}/{n_spots}")

    # Step 2: require at least min_cluster_size neighbors also positive
    adjacency_bin = adjacency.copy()
    adjacency_bin.data = np.ones_like(adjacency_bin.data)
    neighbor_count = adjacency_bin @ candidates.astype(int)
    # A candidate is TLS if it has at least min_cluster_size-1 of its neighbors also candidates
    # (including self, so >= min_cluster_size total)
    tls_label = np.zeros(n_spots, dtype=int)
    tls_label[(candidates) & (neighbor_count >= min_cluster_size)] = 1
    print(f"  TLS-positive spots: {tls_label.sum()}")

    # Step 3: cluster TLS spots into individual structures
    tls_cluster_id = np.zeros(n_spots, dtype=int)
    if tls_label.sum() == 0:
        return tls_label, tls_cluster_id

    # Use connected components on the spatial graph restricted to TLS spots
    tls_indices = np.where(tls_label)[0]
    tls_graph = adjacency[tls_indices][:, tls_indices]
    tls_graph.data = np.ones_like(tls_graph.data)

    # BFS connected components
    visited = np.zeros(len(tls_indices), dtype=bool)
    cluster = 0
    for i in range(len(tls_indices)):
        if visited[i]:
            continue
        cluster += 1
        queue = [i]
        visited[i] = True
        while queue:
            node = queue.pop(0)
            tls_cluster_id[tls_indices[node]] = cluster
            neighbors = tls_graph[node].indices
            for nb in neighbors:
                if not visited[nb]:
                    visited[nb] = True
                    queue.append(nb)

    print(f"  TLS clusters: {cluster}, sizes: {np.bincount(tls_cluster_id[tls_cluster_id>0])}")
    return tls_label, tls_cluster_id


def run_one_sample(h5ad_path: Path, output_dir: Path, z_threshold: float, min_cluster_size: int) -> dict:
    sample_id = h5ad_path.stem
    print(f"\n[sample] {sample_id}")
    adata = ad.read_h5ad(h5ad_path)

    q05 = adata.obsm["c2l_ilc_q05"]
    if hasattr(q05, "values"):
        q05 = q05.values
    cell_types = adata.uns["c2l_ilc_cell_types"]

    # 1. TLS scores
    spot_df = compute_tls_score(q05, cell_types)

    # 2. Spatial graph
    positions = adata.obsm["spatial"]
    adjacency = build_spatial_graph(positions)

    # 3. Label TLS spots and clusters
    tls_label, tls_cluster_id = label_tls_spots(
        spot_df, adjacency, z_threshold=z_threshold, min_cluster_size=min_cluster_size
    )

    # 4. Save results
    sample_out = output_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)

    result = spot_df.copy()
    result["tls_label"] = tls_label
    result["tls_cluster"] = tls_cluster_id
    result.index = adata.obs_names
    result.to_csv(sample_out / "tls_results.csv")

    # Save TLS-positive spot summaries
    tls_spots = result[result["tls_label"] == 1]
    if len(tls_spots) > 0:
        tls_spots.to_csv(sample_out / "tls_positive_spots.csv")

    # Write back to h5ad
    adata.obs["tls_label"] = tls_label
    adata.obs["tls_cluster"] = tls_cluster_id
    adata.obs["tls_zscore"] = result["tls_zscore"].values
    adata.obs["tls_sum"] = result["tls_sum"].values
    adata.write(sample_out / f"{sample_id}_tls.h5ad", compression="gzip")

    return {
        "sample_id": sample_id,
        "n_spots": len(result),
        "tls_positive": int(tls_label.sum()),
        "tls_clusters": int(tls_cluster_id.max()),
    }


def main():
    parser = argparse.ArgumentParser(description="TLS identification from ILC-refined cell2location results")
    parser.add_argument("--anndata-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--z-threshold", type=float, default=1.0)
    parser.add_argument("--min-cluster-size", type=int, default=3)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    h5ad_paths = sorted(args.anndata_dir.glob("*.h5ad"))
    if args.limit > 0:
        h5ad_paths = h5ad_paths[: args.limit]

    rows = []
    for p in h5ad_paths:
        rows.append(run_one_sample(p, args.output, args.z_threshold, args.min_cluster_size))

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output / "tls_summary.csv", index=False)
    n_tls = (summary["tls_positive"] > 0).sum()
    print(f"\n[done] {len(rows)} samples, {n_tls} with TLS-positive spots")


if __name__ == "__main__":
    main()
