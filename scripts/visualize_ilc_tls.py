"""Plot ILC-dominant TLS spots + ILC1/ILC2/ILC3 — each as separate image, then stitch 4x5 grids"""
import argparse
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

COLS, ROWS = 4, 5
N_PER_GRID = 20
TW, TH = 800, 800
ILC_TYPES = ["ILC1", "ILC2", "ILC3"]


def plot_one(sample_dir: Path, h5ad_dir: Path, out_dir: Path) -> None:
    sid = sample_dir.name
    tls_csv = sample_dir / "tls_spot_scores_official_relaxed.csv"
    h5ad_path = h5ad_dir / f"{sid}.h5ad"
    if not tls_csv.exists() or not h5ad_path.exists():
        return

    tls = pd.read_csv(tls_csv)
    adata = ad.read_h5ad(h5ad_path)
    q05 = adata.obsm["c2l_ilc_q05"]
    if hasattr(q05, "values"):
        q05 = q05.values
    ct = list(adata.uns["c2l_ilc_cell_types"])

    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]
    tls_mask = (tls["TLS.region"] == "TLS").values
    dominant = np.array(ct)[q05.argmax(axis=1)]
    # Stricter definition: argmax AND abundance >= sample P75
    p75 = {c: np.percentile(q05[:, ct.index(c)], 75) for c in ILC_TYPES}
    strict_dom = {}
    for c in ILC_TYPES:
        idx = ct.index(c)
        is_top = (dominant == c)
        above_p75 = q05[:, idx] >= p75[c]
        strict_dom[c] = is_top & above_p75
    any_strict = np.any([strict_dom[c] for c in ILC_TYPES], axis=0)

    # --- Panel 1: ILC-dominant spots in TLS (strict: rank=1 AND >= P75) ---
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(x, y, c="lightgray", s=1, rasterized=True)
    colors = {"ILC1": "#e41a1c", "ILC2": "#377eb8", "ILC3": "#4daf4a"}
    has_any = False
    for c in ILC_TYPES:
        mask = tls_mask & strict_dom[c]
        n_total = (tls_mask & (dominant == c)).sum()  # old count for comparison
        if mask.sum() > 0:
            has_any = True
            ax.scatter(x[mask], y[mask], c=colors[c], s=25, edgecolors="black",
                       linewidths=0.5, label=f"{c} ({mask.sum()}/{n_total})", zorder=5)
    if has_any:
        ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"{sid} | ILC-dominant in TLS (rank=1 & >=P75)", fontsize=10)
    ax.axis("off")
    fig.savefig(out_dir / sid / "ilc_tls_dominant.jpg", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Panels 2-4: ILC1, ILC2, ILC3 abundance ---
    for c in ILC_TYPES:
        fig, ax = plt.subplots(figsize=(8, 7))
        idx = ct.index(c)
        vals = q05[:, idx]
        s = ax.scatter(x, y, c=vals, s=2, cmap="YlOrRd", rasterized=True, vmin=0, vmax=np.percentile(vals, 99))
        # TLS boundary
        ax.scatter(x[tls_mask], y[tls_mask], s=20, facecolors="none",
                   edgecolors="red", linewidths=0.8)
        # Strict ILC-dominant TLS spots
        dom_mask = tls_mask & strict_dom[c]
        if dom_mask.sum() > 0:
            ax.scatter(x[dom_mask], y[dom_mask], s=30, marker="*",
                       c="blue", edgecolors="black", linewidths=0.3, zorder=6)
        ax.set_title(f"{sid} | {c}  [P75={p75[c]:.3f}]", fontsize=10)
        ax.axis("off")
        plt.colorbar(s, ax=ax, shrink=0.7)
        fig.savefig(out_dir / sid / f"{c.lower()}_abundance.jpg", dpi=150, bbox_inches="tight")
        plt.close(fig)


def stitch_grids(jpgs: list[Path], out_path: Path, label: str) -> None:
    for g in range(0, len(jpgs), N_PER_GRID):
        batch = jpgs[g:g + N_PER_GRID]
        gi = g // N_PER_GRID + 1
        cells = []
        for f in batch:
            img = Image.open(f)
            img.thumbnail((TW, TH), Image.LANCZOS)
            canvas = Image.new("RGB", (TW, TH), "white")
            canvas.paste(img, ((TW - img.width) // 2, (TH - img.height) // 2))
            cells.append(canvas)
        while len(cells) < N_PER_GRID:
            cells.append(Image.new("RGB", (TW, TH), "white"))
        grid = Image.new("RGB", (COLS * TW, ROWS * TH), "white")
        for i, img in enumerate(cells):
            grid.paste(img, ((i % COLS) * TW, (i // COLS) * TH))
        out = out_path / f"{label}_grid_{gi:02d}.jpg"
        grid.save(out, quality=85)
        print(f"  {label} grid {gi}: {len(batch)} imgs -> {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tls-dir", default=r"E:\GBM\results\tls_official_cut01", type=Path)
    parser.add_argument("--h5ad-dir", default=r"E:\GBM\spatial_data_visium\spatial_data_visium\anndata_with_ilc", type=Path)
    parser.add_argument("--out-dir", default=r"E:\GBM\results\tls_official_cut01", type=Path)
    args = parser.parse_args()

    sample_dirs = sorted(d for d in args.tls_dir.iterdir() if d.is_dir())
    print(f"Plotting {len(sample_dirs)} samples...")
    for sd in sample_dirs:
        plot_one(sd, args.h5ad_dir, args.out_dir)

    features = ["ilc_tls_dominant", "ilc1_abundance", "ilc2_abundance", "ilc3_abundance"]
    for feat in features:
        jpgs = sorted(args.out_dir.glob(f"*/{feat}.jpg"))
        stitch_grids(jpgs, args.out_dir, feat)

    print("Done!")


if __name__ == "__main__":
    main()
