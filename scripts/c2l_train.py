"""Minimal cell2location training + spatial decomp for mgh258"""
import sys, os, warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import numpy as np
import pandas as pd
import cell2location

OUT_DIR = "E:/GBM/results/cell2loc"

# Load reference
print("Loading reference...")
adata_ref = sc.read_h5ad(f"{OUT_DIR}/ref_merged.h5ad")
if adata_ref.raw is not None:
    adata_ref = adata_ref.raw.to_adata()
del adata_ref.uns
sc.pp.filter_genes(adata_ref, min_cells=3)
adata_ref = adata_ref[:, ~adata_ref.var_names.str.startswith('MT-')]
print(f"Reference: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")

# Train
print("\nTraining model...")
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, labels_key='cell_type')
mod = cell2location.models.RegressionModel(adata_ref)
mod.train(max_epochs=150, batch_size=2500, accelerator='cpu')

# Save model
mod.save(f"{OUT_DIR}/ref_model", overwrite=True)
print("Model saved")

# Export signatures
print("Exporting signatures...")
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
adata_ref.write(f"{OUT_DIR}/ref_signatures.h5ad")
print("Signatures saved")

# Load spatial
print("\nLoading spatial...")
adata_sp = sc.read_h5ad(f"{OUT_DIR}/mgh258.h5ad")
if adata_sp.raw is not None:
    adata_sp = adata_sp.raw.to_adata()
del adata_sp.uns
adata_sp = adata_sp[:, ~adata_sp.var_names.str.startswith('MT-')]

# Shared genes
shared = list(set(adata_ref.var_names) & set(adata_sp.var_names))
print(f"Shared genes: {len(shared)}")
adata_ref_sub = adata_ref[:, shared].copy()
adata_sp_sub = adata_sp[:, shared].copy()

# Deconvolve
print("\nRunning cell2location spatial decomposition...")
adata_sp_sub = cell2location.run_cell2location(
    adata_ref=adata_ref_sub,
    adata_sp=adata_sp_sub,
    N_cells_per_location=8,
    detection_alpha=20,
    max_epochs=30000,
    batch_size=None
)

# Save
adata_sp_sub.write(f"{OUT_DIR}/mgh258_cell2loc.h5ad")
df = pd.DataFrame(adata_sp_sub.obsm['q05_cell_abundance_w_sf'],
                  index=adata_sp_sub.obs_names,
                  columns=adata_sp_sub.uns['mod']['factor_names'])
df.to_csv(f"{OUT_DIR}/mgh258_cell_proportions.csv")
print(f"\nDone! Proportions: {df.shape}")
print(df.head())
