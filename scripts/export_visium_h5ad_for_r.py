from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from scipy import io as spio
from scipy import sparse


def write_tsv_gz(lines: list[str], path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")


def write_csv_gz(df: pd.DataFrame, path: Path, index: bool = True) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        df.to_csv(handle, index=index)


def write_mtx_gz(matrix: sparse.spmatrix, path: Path) -> None:
    with gzip.open(path, "wb") as handle:
        spio.mmwrite(handle, matrix)


def sanitize_image(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image)
    if arr.dtype != np.uint8:
        if np.issubdtype(arr.dtype, np.floating):
            scale = 255.0 if arr.max(initial=0) <= 1.0 else 1.0
            arr = np.clip(arr * scale, 0, 255).astype(np.uint8)
        else:
            arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def get_spatial_block(adata: ad.AnnData) -> tuple[str, dict]:
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, dict) or not spatial:
        raise ValueError("Missing adata.uns['spatial']")
    library_id = next(iter(spatial))
    return library_id, spatial[library_id]


def export_one(h5ad_path: Path, out_root: Path) -> dict:
    sample_id = h5ad_path.stem
    sample_dir = out_root / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)

    obs = adata.obs.copy()
    obs.index.name = "barcode"
    var = adata.var.copy()
    var.index.name = "feature_name"

    library_id, spatial_block = get_spatial_block(adata)
    images = spatial_block.get("images", {})
    scalefactors = spatial_block.get("scalefactors", {})
    metadata = spatial_block.get("metadata", {})

    gene_mask = np.asarray(var["feature_types"] == "Gene Expression")
    gene_var = var.loc[gene_mask].copy()
    gene_x = x[:, gene_mask]

    gene_dir = sample_dir / "gene_expression"
    all_dir = sample_dir / "all_features"
    spatial_dir = sample_dir / "spatial"
    gene_dir.mkdir(exist_ok=True)
    all_dir.mkdir(exist_ok=True)
    spatial_dir.mkdir(exist_ok=True)

    barcodes = obs.index.astype(str).tolist()
    write_tsv_gz(barcodes, gene_dir / "barcodes.tsv.gz")
    write_tsv_gz(barcodes, all_dir / "barcodes.tsv.gz")

    gene_features = [
        f"{feature}\t{feature}\tGene Expression"
        for feature in gene_var.index.astype(str).tolist()
    ]
    write_tsv_gz(gene_features, gene_dir / "features.tsv.gz")

    all_features = [
        f"{idx}\t{idx}\t{ftype}"
        for idx, ftype in zip(var.index.astype(str), var["feature_types"].astype(str))
    ]
    write_tsv_gz(all_features, all_dir / "features.tsv.gz")

    # Matrix Market for 10x/Seurat convention: features x barcodes.
    write_mtx_gz(gene_x.T.tocoo(), gene_dir / "matrix.mtx.gz")
    write_mtx_gz(x.T.tocoo(), all_dir / "matrix.mtx.gz")

    write_csv_gz(obs, sample_dir / "obs.csv.gz", index=True)
    write_csv_gz(var, sample_dir / "var.csv.gz", index=True)
    write_csv_gz(gene_var, sample_dir / "gene_expression_var.csv.gz", index=True)

    spatial_coords = pd.DataFrame(
        adata.obsm["spatial"],
        index=obs.index,
        columns=["pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    spatial_coords.insert(0, "array_col", obs.get("array_col"))
    spatial_coords.insert(0, "array_row", obs.get("array_row"))
    spatial_coords.insert(0, "in_tissue", obs.get("in_tissue"))
    spatial_coords.insert(0, "barcode", obs.index.astype(str))
    spatial_coords.to_csv(spatial_dir / "tissue_positions.csv", index=False)
    write_csv_gz(spatial_coords.set_index("barcode"), spatial_dir / "coordinates.csv.gz", index=True)

    with open(spatial_dir / "scalefactors_json.json", "w", encoding="utf-8") as handle:
        json.dump(scalefactors, handle, ensure_ascii=False, indent=2)

    with open(spatial_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)

    for image_name in ("hires", "lowres"):
        if image_name in images:
            Image.fromarray(sanitize_image(images[image_name])).save(
                spatial_dir / f"tissue_{image_name}_image.png"
            )

    with open(sample_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "sample_id": sample_id,
                "library_id": library_id,
                "source_h5ad": str(h5ad_path),
                "n_spots": int(adata.n_obs),
                "n_all_features": int(adata.n_vars),
                "n_gene_expression_features": int(gene_mask.sum()),
                "feature_types": var["feature_types"].value_counts().to_dict(),
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )

    return {
        "sample_id": sample_id,
        "n_spots": int(adata.n_obs),
        "n_all_features": int(adata.n_vars),
        "n_gene_features": int(gene_mask.sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export Visium h5ad files into R-reconstructable directories."
    )
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--pattern",
        default="*.h5ad",
        help="Glob pattern for h5ad files inside input_dir.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_files = sorted(input_dir.glob(args.pattern))
    if not h5ad_files:
        raise SystemExit(f"No files matched {args.pattern} under {input_dir}")

    summary = []
    for h5ad_path in h5ad_files:
        print(f"[export] {h5ad_path.name}")
        summary.append(export_one(h5ad_path, output_dir))

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / "export_summary.csv", index=False)
    print(f"[done] exported {len(summary_df)} files to {output_dir}")


if __name__ == "__main__":
    main()
