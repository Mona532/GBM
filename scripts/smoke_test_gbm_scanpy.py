from __future__ import annotations

import ast
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


ROOT = Path(r"E:\GBM")
DATA_DIR = ROOT / "GBM_DATA" / "5sample_final"
RESULTS_DIR = ROOT / "results"

NOTEBOOK_PY = ROOT / "GBM_scanpy.py"
INPUT_H5AD = DATA_DIR / "GBM_adata_6_11.h5ad"
ANNO_H5AD = DATA_DIR / "GBM_5sample_anno.h5ad"
ILC_H5AD = DATA_DIR / "GBM_ilc.h5ad"
MARKERS_CSV = DATA_DIR / "markers_df.csv"
SMOKE_PNG = RESULTS_DIR / "gbm_scanpy_smoke_umap.png"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_notebook_script() -> None:
    require(NOTEBOOK_PY.exists(), f"Missing notebook script: {NOTEBOOK_PY}")
    source = NOTEBOOK_PY.read_text(encoding="utf-8")
    ast.parse(source, filename=str(NOTEBOOK_PY))
    print(f"[ok] Parsed notebook script: {NOTEBOOK_PY.name}")


def check_input_h5ad() -> None:
    require(INPUT_H5AD.exists(), f"Missing input h5ad: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    require(adata.n_obs > 1000, "Input h5ad has too few cells")
    require(adata.n_vars > 1000, "Input h5ad has too few genes")
    for col in ["sample", "Source", "leiden_r0.5", "leiden_r0.6"]:
        require(col in adata.obs.columns, f"Input h5ad missing obs column: {col}")
    for key in ["X_pca", "X_pca_harmony", "X_umap"]:
        require(key in adata.obsm, f"Input h5ad missing embedding: {key}")
    print(f"[ok] Input h5ad: {adata.n_obs:,} cells x {adata.n_vars:,} genes")


def check_markers_csv() -> None:
    require(MARKERS_CSV.exists(), f"Missing markers CSV: {MARKERS_CSV}")
    markers = pd.read_csv(MARKERS_CSV)
    required_cols = {
        "cluster",
        "names",
        "logfoldchanges",
        "pvals_adj",
        "pct_nz_group",
        "pct_nz_reference",
    }
    require(required_cols.issubset(markers.columns), "Markers CSV is missing required columns")
    require(not markers.empty, "Markers CSV is empty")
    require(markers["cluster"].nunique() >= 3, "Markers CSV has too few clusters")
    print(f"[ok] Markers CSV: {len(markers):,} rows across {markers['cluster'].nunique()} clusters")


def check_annotated_outputs() -> None:
    require(ANNO_H5AD.exists(), f"Missing annotated h5ad: {ANNO_H5AD}")
    require(ILC_H5AD.exists(), f"Missing ILC h5ad: {ILC_H5AD}")

    anno = sc.read_h5ad(ANNO_H5AD)
    ilc = sc.read_h5ad(ILC_H5AD)

    for col in ["sample", "broad_anno", "leiden_r0.5", "leiden_r0.6"]:
        require(col in anno.obs.columns, f"Annotated h5ad missing obs column: {col}")
    require("X_umap" in anno.obsm, "Annotated h5ad missing UMAP embedding")

    for col in ["sample", "leiden_r1.0", "ilc_subtype"]:
        require(col in ilc.obs.columns, f"ILC h5ad missing obs column: {col}")
    for key in ["X_umap", "X_pca_harmony_ilc"]:
        require(key in ilc.obsm, f"ILC h5ad missing embedding: {key}")
    require("lognorm" in ilc.layers, "ILC h5ad missing lognorm layer")
    require(ilc.obs["ilc_subtype"].nunique() >= 2, "ILC h5ad has too few ILC subtypes")

    print(f"[ok] Annotated h5ad: {anno.n_obs:,} cells")
    print(f"[ok] ILC h5ad: {ilc.n_obs:,} cells across {ilc.obs['ilc_subtype'].nunique()} subtypes")


def write_smoke_plot() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ilc = sc.read_h5ad(ILC_H5AD)
    fig = sc.pl.umap(
        ilc,
        color="ilc_subtype",
        show=False,
        return_fig=True,
        title="GBM ILC Smoke Test",
    )
    fig.savefig(SMOKE_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    require(SMOKE_PNG.exists(), f"Failed to write smoke plot: {SMOKE_PNG}")
    print(f"[ok] Wrote smoke plot: {SMOKE_PNG}")


def main() -> None:
    check_notebook_script()
    check_input_h5ad()
    check_markers_csv()
    check_annotated_outputs()
    write_smoke_plot()
    print("[done] GBM scanpy smoke test passed")


if __name__ == "__main__":
    main()
