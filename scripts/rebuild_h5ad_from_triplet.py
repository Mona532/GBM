"""
Rebuild h5ad from triplet files (counts.mtx + barcodes.tsv + features.tsv + metadata.csv)
- X = raw counts (sparse, cells × genes)
- obs = cell metadata from metadata.csv
- var = gene names from features.tsv
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix, issparse

TRIPLET_DIR = "E:/GBM/GBM_DATA/5sample_integration/triplet"
OUT_H5AD     = "E:/GBM/GBM_DATA/5sample_integration/GBM_5samples_filtered_clean.h5ad"

# ── 1. Read count matrix (R format: genes × cells, MM format as written) ──
print("Loading counts.mtx...")
counts = mmread(f"{TRIPLET_DIR}/counts.mtx")  # genes × cells (dgTMatrix)
print(f"  Raw shape (genes × cells): {counts.shape}")

# Transpose to cells × genes for anndata
counts = counts.T.tocsr()  # 18488 cells × 25622 genes
print(f"  Transposed (cells × genes): {counts.shape}")
print(f"  Sparsity: {100 * (1 - counts.nnz / (counts.shape[0] * counts.shape[1])):.2f}%")

# ── 2. Read barcodes ──
barcodes = pd.read_csv(f"{TRIPLET_DIR}/barcodes.tsv", header=None, names=["barcode"])
print(f"  Barcodes: {len(barcodes)}")

# ── 3. Read features ──
features = pd.read_csv(f"{TRIPLET_DIR}/features.tsv", sep="\t", header=None,
                       names=["gene_id", "gene_name", "feature_type"])
print(f"  Features: {len(features)}")

# ── 4. Read metadata ──
meta = pd.read_csv(f"{TRIPLET_DIR}/metadata.csv", index_col=0)
print(f"  Metadata: {meta.shape[1]} cols × {meta.shape[0]} rows")

# ── 5. Build AnnData ──
adata = sc.AnnData(
    X=counts,
    obs=meta,
    var=features.copy()
)
adata.obs_names = barcodes["barcode"].values
adata.var_names = features["gene_name"].values
adata.var_names_make_unique()

print(f"\n✓ AnnData built: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"  X range: {adata.X.min():.1f} – {adata.X.max():.1f}")
print(f"  X mean: {adata.X.mean():.2f}")

# ── 6. Save ──
adata.write(OUT_H5AD, compression="gzip")
print(f"\n✓ Saved to: {OUT_H5AD}")
print(f"  File size: {__import__('os').path.getsize(OUT_H5AD) / 1e6:.1f} MB")
