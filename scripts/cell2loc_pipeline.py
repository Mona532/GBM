"""
Cell2location deconvolution for GBM Visium data (merged cell types)
SpaLinker paper workflow: scRNA-seq reference → Cell2location → spot proportions → SpaLinker TLS
"""
import sys, os, warnings
import scanpy as sc
import numpy as np
import pandas as pd
import cell2location
import torch
warnings.filterwarnings("ignore", category=UserWarning)

OUT_DIR = "E:/GBM/results/cell2loc"

# ===================================================================
# Step 1: Load merged reference
# ===================================================================
print("=" * 60)
print("Loading merged reference (12 broad cell types)")
print("=" * 60)
adata_ref = sc.read_h5ad(f"{OUT_DIR}/ref_merged.h5ad")
print(f"Reference: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")

# Get raw counts (cell2location needs them)
if adata_ref.raw is not None:
    adata_ref = adata_ref.raw.to_adata()
del adata_ref.uns  # clean up

# Filter: keep only genes expressed in >=3 cells
sc.pp.filter_genes(adata_ref, min_cells=3)

# Remove mitochondrial genes from reference
adata_ref = adata_ref[:, ~adata_ref.var_names.str.startswith('MT-')]

print(f"After filtering: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")
print(f"Cell types ({adata_ref.obs['cell_type'].nunique()}):")
print(adata_ref.obs['cell_type'].value_counts())

# ===================================================================
# Step 2: Train reference signatures
# ===================================================================
print("\n" + "=" * 60)
print("Training RegressionModel (GPU)")
print("=" * 60)

cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    labels_key='cell_type'
)

mod = cell2location.models.RegressionModel(adata_ref)
mod.train(max_epochs=150, batch_size=2500, accelerator='gpu')

mod.save(f"{OUT_DIR}/ref_model", overwrite=True)

# Export posterior
adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)
adata_ref.write(f"{OUT_DIR}/ref_signatures.h5ad")
print("Reference signatures saved")

# ===================================================================
# Step 3: Load spatial data
# ===================================================================
print("\n" + "=" * 60)
print("Loading spatial Visium data")
print("=" * 60)
adata_sp = sc.read_h5ad(f"{OUT_DIR}/mgh258.h5ad")
print(f"Spatial: {adata_sp.n_obs} spots x {adata_sp.n_vars} genes")

# Use raw counts for spatial too
if adata_sp.raw is not None:
    adata_sp = adata_sp.raw.to_adata()
del adata_sp.uns

# Filter mitochondrial genes
adata_sp = adata_sp[:, ~adata_sp.var_names.str.startswith('MT-')]

# Find shared genes
shared = list(set(adata_ref.var_names) & set(adata_sp.var_names))
print(f"Shared genes: {len(shared)}")
adata_ref_sub = adata_ref[:, shared].copy()
adata_sp_sub  = adata_sp[:, shared].copy()

# ===================================================================
# Step 4: Run cell2location
# ===================================================================
print("\n" + "=" * 60)
print("Running cell2location spatial decomposition (GPU)")
print("=" * 60)

adata_sp_sub = cell2location.run_cell2location(
    adata_ref=adata_ref_sub,
    adata_sp=adata_sp_sub,
    N_cells_per_location=8,
    detection_alpha=20,
    max_epochs=30000,
    batch_size=None
)

adata_sp_sub.write(f"{OUT_DIR}/mgh258_cell2loc.h5ad")

# ===================================================================
# Step 5: Export cell proportions for SpaLinker
# ===================================================================
cell_abundances = adata_sp_sub.obsm['q05_cell_abundance_w_sf']
cell_types = adata_sp_sub.uns['mod']['factor_names']
print(f"\nCell type proportions shape: {cell_abundances.shape}")
print(f"Cell types: {list(cell_types)}")

df = pd.DataFrame(cell_abundances, index=adata_sp_sub.obs_names, columns=cell_types)
# Rename columns to match SpaLinker naming: paper uses "Plasma/B cells" but we have "B cell"
# "B cell" and "T cell" are the key names CalCellCodis will use
df.to_csv(f"{OUT_DIR}/mgh258_cell_proportions.csv")
print(f"\nSaved to: {OUT_DIR}/mgh258_cell_proportions.csv")
print("Done!")
