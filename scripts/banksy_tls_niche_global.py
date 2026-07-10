"""Global-consistent BANKSY niche analysis for GBM TLS.

The previous niche workflow clustered each sample independently and then
aggregated cluster IDs across samples. KMeans labels are arbitrary per fit,
so that aggregation is not stable. This script fixes that by:

1. Building a BANKSY embedding per sample from compositional c2l abundance.
2. Concatenating all sample embeddings into one pooled matrix.
3. Fitting one global clustering model to assign consistent niche labels.
4. Summarizing TLS enrichment and niche composition with sample-level support.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from banksy.embed_banksy import generate_banksy_matrix
from banksy.initialize_banksy import initialize_banksy
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler

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
K = 5
LAMBDA = 0.2
N_NEIGH = 15
MIN_SHARED = 100
MIN_TLS = 5
MIN_NICHE_SPOTS_PER_SAMPLE = 20


@dataclass
class SampleBundle:
    sample: str
    prop: np.ndarray
    coords: np.ndarray
    tls_mask: np.ndarray
    ct_names: list[str]
    embedding: np.ndarray | None = None
    labels: np.ndarray | None = None


def _to_numpy(x) -> np.ndarray:
    if hasattr(x, "values"):
        x = x.values
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x)


def load_samples() -> list[SampleBundle]:
    bundles: list[SampleBundle] = []

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

        q05 = _to_numpy(adata.obsm["c2l_ilc_q05"]).astype(np.float64)
        row_sum = q05.sum(axis=1, keepdims=True) + 1e-8
        prop = q05 / row_sum
        coords = _to_numpy(adata.obsm["spatial"]).astype(np.float64)
        ct_names = list(adata.uns["c2l_ilc_cell_types"])

        bundles.append(
            SampleBundle(
                sample=h5_path.stem,
                prop=prop,
                coords=coords,
                tls_mask=tls_mask,
                ct_names=ct_names,
            )
        )

    return bundles


def build_banksy_embeddings(bundles: list[SampleBundle]) -> None:
    for i, bundle in enumerate(bundles, start=1):
        ba = ad.AnnData(
            X=bundle.prop.astype(np.float64),
            obs=pd.DataFrame({"sample": [bundle.sample] * bundle.prop.shape[0]}),
            obsm={"spatial": bundle.coords, "spatial_coords": bundle.coords},
        )
        ba.var_names = bundle.ct_names

        bd = initialize_banksy(
            ba,
            coord_keys=("spatial", "spatial", "spatial_coords"),
            num_neighbours=N_NEIGH,
            max_m=0,
            plt_edge_hist=False,
            plt_nbr_weights=False,
            plt_theta=False,
        )
        _, bm = generate_banksy_matrix(ba, bd, [LAMBDA], max_m=0)
        bundle.embedding = _to_numpy(bm.X).astype(np.float32)

        if i % 30 == 0:
            print(f"  BANKSY embeddings: {i}/{len(bundles)}")


def assign_global_niches(bundles: list[SampleBundle]) -> tuple[np.ndarray, StandardScaler]:
    stacked = np.vstack([b.embedding for b in bundles if b.embedding is not None])
    scaler = StandardScaler()
    stacked_scaled = scaler.fit_transform(stacked)

    model = MiniBatchKMeans(
        n_clusters=K,
        random_state=42,
        batch_size=8192,
        n_init=20,
        max_iter=200,
    )
    global_labels = model.fit_predict(stacked_scaled)

    offset = 0
    for bundle in bundles:
        n = bundle.embedding.shape[0]
        bundle.labels = global_labels[offset : offset + n]
        offset += n

    return global_labels, scaler


def summarize_niches(bundles: list[SampleBundle]) -> pd.DataFrame:
    ct_names = bundles[0].ct_names
    all_prop = np.vstack([b.prop for b in bundles])
    all_labels = np.concatenate([b.labels for b in bundles])
    all_tls = np.concatenate([b.tls_mask for b in bundles])

    global_mean = all_prop.mean(axis=0) + 1e-8
    rows = []

    for niche in range(K):
        niche_mask = all_labels == niche
        niche_prop = all_prop[niche_mask]
        niche_mean = niche_prop.mean(axis=0)
        niche_enrich = np.log2((niche_mean + 1e-8) / global_mean)

        tls_count = int((all_tls & niche_mask).sum())
        all_count = int(niche_mask.sum())
        tls_pct = 100.0 * tls_count / max(int(all_tls.sum()), 1)
        all_pct = 100.0 * all_count / max(len(all_labels), 1)
        tls_log2fc = np.log2(((tls_count + 1) / (all_count + 1)) * (len(all_labels) / max(int(all_tls.sum()), 1)))

        n_samples = 0
        for bundle in bundles:
            if int((bundle.labels == niche).sum()) >= MIN_NICHE_SPOTS_PER_SAMPLE:
                n_samples += 1

        row = {ct: niche_enrich[i] for i, ct in enumerate(ct_names)}
        row["niche"] = niche
        row["tls_log2fc"] = tls_log2fc
        row["tls_pct"] = tls_pct
        row["all_pct"] = all_pct
        row["n_spots"] = all_count
        row["n_tls_spots"] = tls_count
        row["n_samples"] = n_samples
        rows.append(row)

    return pd.DataFrame(rows)


def save_spot_assignments(bundles: list[SampleBundle]) -> None:
    assign_dir = OUT / "banksy_global_assignments"
    assign_dir.mkdir(parents=True, exist_ok=True)

    for bundle in bundles:
        df = pd.DataFrame(
            {
                "barcode": np.arange(bundle.labels.shape[0]),
                "niche": bundle.labels,
                "tls": bundle.tls_mask.astype(int),
                "x": bundle.coords[:, 0],
                "y": bundle.coords[:, 1],
            }
        )
        df.to_csv(assign_dir / f"{bundle.sample}.csv", index=False)


def plot_summary(df: pd.DataFrame, ct_names: list[str]) -> None:
    enrich = df[ct_names].to_numpy()
    tls_log2fc = df["tls_log2fc"].to_numpy()

    fig = plt.figure(figsize=(12, 6))

    ax1 = fig.add_axes([0.05, 0.12, 0.38, 0.80])
    vmax = max(abs(enrich.min()), abs(enrich.max()), 0.5)
    im1 = ax1.imshow(enrich.T, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax1.set_xticks(range(K))
    ax1.set_xticklabels(
        [f"N{k}\nn={int(df.loc[k, 'n_spots'])}\nTLS {tls_log2fc[k]:+.2f}" for k in range(K)],
        fontsize=6,
    )
    ax1.set_yticks(range(len(ct_names)))
    ax1.set_yticklabels(ct_names, fontsize=6)
    ax1.set_title(
        f"a  Global BANKSY niches (K={K}, lambda={LAMBDA})",
        fontsize=9,
        fontweight="bold",
        loc="left",
    )
    plt.colorbar(im1, ax=ax1, shrink=0.7).set_label("log2 enrichment", fontsize=6)

    ax2 = fig.add_axes([0.50, 0.12, 0.20, 0.80])
    colors = ["#d7191c" if v > 0 else "#2b83ba" for v in tls_log2fc]
    ax2.bar(range(K), tls_log2fc, color=colors, alpha=0.75, edgecolor="black", linewidth=0.3)
    ax2.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax2.set_xticks(range(K))
    ax2.set_xticklabels([f"N{k}" for k in range(K)], fontsize=7)
    ax2.set_ylabel("TLS enrichment (log2 FC)", fontsize=7)
    ax2.set_title("b  TLS enrichment", fontsize=9, fontweight="bold", loc="left")

    ax3 = fig.add_axes([0.76, 0.12, 0.22, 0.80])
    support = df[["n_samples", "n_tls_spots", "all_pct"]].to_numpy(dtype=float)
    support_z = (support - support.mean(axis=0)) / (support.std(axis=0) + 1e-8)
    im3 = ax3.imshow(support_z.T, aspect="auto", cmap="YlGnBu", vmin=-1.5, vmax=1.5)
    ax3.set_xticks(range(K))
    ax3.set_xticklabels([f"N{k}" for k in range(K)], fontsize=7)
    ax3.set_yticks(range(3))
    ax3.set_yticklabels(["samples", "TLS spots", "all %"], fontsize=6)
    ax3.set_title("c  Support metrics", fontsize=9, fontweight="bold", loc="left")
    plt.colorbar(im3, ax=ax3, shrink=0.7).set_label("z-score", fontsize=6)

    fig.savefig(OUT / "fig_banksy_tls_niche_global.jpg", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bundles = load_samples()
    if not bundles:
        raise RuntimeError("No eligible samples found for BANKSY niche analysis.")

    print(f"Samples: {len(bundles)}")
    build_banksy_embeddings(bundles)
    assign_global_niches(bundles)

    summary = summarize_niches(bundles)
    summary.to_csv(OUT / "banksy_niche_tls_global.csv", index=False)
    save_spot_assignments(bundles)
    plot_summary(summary, bundles[0].ct_names)

    print(f"Saved: {OUT / 'banksy_niche_tls_global.csv'}")
    print(f"Saved: {OUT / 'fig_banksy_tls_niche_global.jpg'}")
    for niche in range(K):
        top_ct = summary.loc[niche, bundles[0].ct_names].idxmax()
        print(
            f"N{niche}: top={top_ct}, TLS_log2FC={summary.loc[niche, 'tls_log2fc']:+.2f}, "
            f"samples={int(summary.loc[niche, 'n_samples'])}"
        )


if __name__ == "__main__":
    main()
