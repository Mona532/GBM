"""Inspect spatial anndata structure"""
import scanpy as sc
import pandas as pd
import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

DATA = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"

# Read README with error handling
try:
    with open(f"{DATA}/README.md", encoding='utf-8', errors='replace') as f:
        print(f.read()[:2000])
except Exception as e:
    print(f"README error: {e}")

# Read one small sample
fp = f"{DATA}/AT3-BRA5-FO-1_1.h5ad"
print(f"\n{'='*60}")
print(f"Reading: {fp}")
adata = sc.read_h5ad(fp)
print(f"Shape: {adata.shape}")
print(f"\nobs columns ({len(adata.obs.columns)}):")
for c in adata.obs.columns:
    n = adata.obs[c].nunique()
    dtype = adata.obs[c].dtype
    if n < 50 and n > 1:
        print(f"  [{c}] dtype={dtype}, n={n}")
        vc = adata.obs[c].value_counts()
        for k, v in vc.items():
            print(f"    {k}: {v}")
    else:
        print(f"  [{c}] dtype={dtype}, n={n} (skip detail)")

print(f"\nobsm keys: {list(adata.obsm.keys())}")
if 'spatial' in adata.obsm:
    print(f"  spatial shape: {adata.obsm['spatial'].shape}")
    print(f"  spatial[:3]: {adata.obsm['spatial'][:3]}")

print(f"\nuns keys: {list(adata.uns.keys())}")

# Check X type
print(f"\nX type: {type(adata.X)}")
if hasattr(adata.X, 'shape'):
    print(f"X shape: {adata.X.shape}")
