"""Extract cell2location proportions from anndata and check cell type names"""
import scanpy as sc
import pandas as pd
import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

DATA = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"
fp = f"{DATA}/AT3-BRA5-FO-1_1.h5ad"

adata = sc.read_h5ad(fp)

# Check feature types
print("=== Feature types ===")
print(adata.var['feature_type'].value_counts())

# Extract cell state abundances
cell_abund = adata[:, adata.var['feature_type'] == 'Cell state abundances']
print(f"\nCell state abundance features: {cell_abund.shape[1]}")

# Show all cell state names
cell_names = cell_abund.var_names.tolist()
print(f"\nAll {len(cell_names)} cell state names:")
for i, name in enumerate(cell_names):
    print(f"  [{i}] {name}")

# Convert to dense dataframe
abund_dense = pd.DataFrame(
    cell_abund.X.toarray() if hasattr(cell_abund.X, 'toarray') else cell_abund.X,
    index=cell_abund.obs_names,
    columns=cell_names
)
print(f"\nAbundance matrix: {abund_dense.shape}")
print(f"Range: [{abund_dense.values.min():.4f}, {abund_dense.values.max():.4f}]")
print(f"Head:")
print(abund_dense.head())

# Check spatial coordinates
print(f"\nSpatial coords shape: {adata.obsm['spatial'].shape}")
print(f"Sample: {adata.obs['sample_name'].iloc[0]}")
