from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


TLS_CELL_KEYWORDS = [
    "B cells",
    "Plasma",
    "CD4",
    "CD8",
    "T cells",
    "T cell",
    "Dendritic",
    "lymphocyte",
]

TLS_NICHE_KEYWORDS = [
    "Immune",
    "Vasculature",
]

TLS_HISTO_KEYWORDS = [
    "Cellular tumor",
    "Leading edge",
    "Necrosis",
    "Perinecrotic",
    "Hyperplastic blood vessels",
    "Microvascular proliferation",
]


def pick_matches(values: list[str], keywords: list[str]) -> list[str]:
    matched: list[str] = []
    for value in values:
        lower_value = value.lower()
        if any(keyword.lower() in lower_value for keyword in keywords):
            matched.append(value)
    return matched


def scan_one(path: Path) -> dict:
    adata = ad.read_h5ad(path)
    var = adata.var.copy()
    feature_types = var["feature_types"].astype(str)

    gene_features = var.index[feature_types == "Gene Expression"].astype(str).tolist()
    cell_states = var.index[feature_types == "Cell state abundances"].astype(str).tolist()
    niches = var.index[feature_types == "Spatial niche abundances"].astype(str).tolist()
    histo = var.index[feature_types == "Histopath annotation overlap"].astype(str).tolist()

    tls_cells = pick_matches(cell_states, TLS_CELL_KEYWORDS)
    tls_niches = pick_matches(niches, TLS_NICHE_KEYWORDS)
    tls_histo = pick_matches(histo, TLS_HISTO_KEYWORDS)

    return {
        "sample_id": path.stem,
        "n_spots": adata.n_obs,
        "n_features": adata.n_vars,
        "n_gene_expression": len(gene_features),
        "n_cell_state": len(cell_states),
        "n_spatial_niche": len(niches),
        "n_histopath": len(histo),
        "has_spatial_coords": "spatial" in adata.obsm,
        "has_spatial_uns": "spatial" in adata.uns,
        "has_b_cells": "B cells" in cell_states,
        "has_plasma_cells": any("plasma" in x.lower() for x in cell_states),
        "has_cd4_t": any("cd4" in x.lower() for x in cell_states),
        "has_cd8_t": any("cd8" in x.lower() for x in cell_states),
        "has_any_t_cells": any("t cell" in x.lower() or "t cells" in x.lower() for x in cell_states),
        "has_dendritic": any("dendritic" in x.lower() for x in cell_states),
        "has_ambiguous_lymphocyte": any("lymphocyte" in x.lower() for x in cell_states),
        "tls_cell_features": "; ".join(tls_cells),
        "tls_niche_features": "; ".join(tls_niches),
        "tls_histopath_features": "; ".join(tls_histo),
        "all_cell_state_features": "; ".join(cell_states),
        "all_niche_features": "; ".join(niches),
        "all_histopath_features": "; ".join(histo),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Scan Visium h5ad files for TLS-relevant features.")
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("*.h5ad"))
    if not files:
        raise SystemExit(f"No .h5ad files found under {args.input_dir}")

    rows = []
    for path in files:
        print(f"[scan] {path.name}")
        rows.append(scan_one(path))

    df = pd.DataFrame(rows).sort_values("sample_id")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)

    quick = {
        "n_samples": int(df.shape[0]),
        "samples_with_b_cells": int(df["has_b_cells"].sum()),
        "samples_with_plasma_cells": int(df["has_plasma_cells"].sum()),
        "samples_with_cd4_t": int(df["has_cd4_t"].sum()),
        "samples_with_cd8_t": int(df["has_cd8_t"].sum()),
        "samples_with_dendritic": int(df["has_dendritic"].sum()),
        "samples_with_histopath": int((df["n_histopath"] > 0).sum()),
    }
    print(pd.Series(quick).to_string())


if __name__ == "__main__":
    main()
