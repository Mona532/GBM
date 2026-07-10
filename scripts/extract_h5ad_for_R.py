"""Extract h5ad data to R-friendly format: mtx + barcodes + features + coords + cell_props"""
import scanpy as sc
import scipy.sparse as sp
import scipy.io
import pandas as pd
import numpy as np
import os, sys

DATA_DIR = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"
OUT_DIR  = "E:/GBM/results/spalinker_input"
os.makedirs(OUT_DIR, exist_ok=True)

# Process first sample
fp = os.path.join(DATA_DIR, "AT3-BRA5-FO-1_1.h5ad")
print(f"Reading {fp}...")
adata = sc.read_h5ad(fp)
print(f"Shape: {adata.shape}")
print(f"obs columns: {list(adata.obs.columns)}")
print(f"var columns: {list(adata.var.columns)}")
print(f"obsm keys: {list(adata.obsm.keys())}")

# Extract gene expression
gene_mask = adata.var['feature_type'] == 'Gene Expression'
print(f"Gene Expression features: {gene_mask.sum()}")
adata_gene = adata[:, gene_mask].copy()

# Extract cell state abundances
cell_mask = adata.var['feature_type'] == 'Cell state abundances'
print(f"Cell state abundance features: {cell_mask.sum()}")
adata_cell = adata[:, cell_mask].copy()

# Get cell state names
cell_names = list(adata_cell.var_names)
print("Cell state names:")
for i, n in enumerate(cell_names):
    print(f"  [{i}] {n}")

# Save gene expression as mtx
print("\nSaving gene expression...")
gene_counts = adata_gene.X
if sp.issparse(gene_counts):
    sp.io.mmwrite(f"{OUT_DIR}/gene_counts.mtx", gene_counts)
else:
    sp.io.mmwrite(f"{OUT_DIR}/gene_counts.mtx", sp.csr_matrix(gene_counts))

# Save barcodes
with open(f"{OUT_DIR}/barcodes.tsv", 'w') as f:
    for b in adata.obs_names:
        f.write(f"{b}\n")

# Save gene names
with open(f"{OUT_DIR}/genes.tsv", 'w') as f:
    for g in adata_gene.var_names:
        f.write(f"{g}\n")

# Save cell state names
with open(f"{OUT_DIR}/cell_state_names.tsv", 'w') as f:
    for c in cell_names:
        f.write(f"{c}\n")

# Save spatial coordinates
spatial = adata.obsm['spatial']
np.savetxt(f"{OUT_DIR}/spatial_coords.tsv", spatial, delimiter='\t', header='x\ty', comments='')

# Save cell state abundance matrix (dense)
cell_dense = adata_cell.X.toarray() if sp.issparse(adata_cell.X) else adata_cell.X
np.savetxt(f"{OUT_DIR}/cell_props.tsv", cell_dense, delimiter='\t')

# Also check B/T cell state identification
b_states = [n for n in cell_names if any(k in n.lower() for k in ['b cell','b_cell','plasma','bcell','b lymph'])]
t_states = [n for n in cell_names if any(k in n.lower() for k in ['t cell','t_cell','cd4','cd8','tcell','t lymph'])]
print(f"\nB/Plasma states: {b_states}")
print(f"T cell states: {t_states}")

print(f"\nAll data saved to {OUT_DIR}")
