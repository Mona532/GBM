"""Simple extraction: gene expression (mtx) + cell abundances (tsv) + spatial coords"""
import scanpy as sc
import scipy.sparse as sp
import scipy.io
import numpy as np
import os, glob

DATA = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"
OUT  = "E:/GBM/results/spalinker_input"
os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(f"{DATA}/*.h5ad"))
print(f"Files: {len(files)}")

for i, fp in enumerate(files):
    tag = os.path.basename(fp).replace('.h5ad', '')
    odir = os.path.join(OUT, tag)
    os.makedirs(odir, exist_ok=True)
    if os.path.exists(os.path.join(odir, 'done.txt')):
        continue

    try:
        adata = sc.read_h5ad(fp)
        # Gene Expression
        gmask = adata.var['feature_types'] == 'Gene Expression'
        sp.io.mmwrite(os.path.join(odir, 'counts.mtx'), adata[:, gmask].X)
        np.savetxt(os.path.join(odir, 'genes.txt'), adata.var_names[gmask].values, fmt='%s')
        # Cell state abundances
        cmask = adata.var['feature_types'] == 'Cell state abundances'
        if cmask.sum() > 0:
            cdat = adata[:, cmask].X
            if sp.issparse(cdat): cdat = cdat.toarray()
            np.savetxt(os.path.join(odir, 'cell_abund.tsv'), cdat, delimiter='\t')
            np.savetxt(os.path.join(odir, 'cell_names.txt'), adata.var_names[cmask].values, fmt='%s')
        # Barcodes + coords
        with open(os.path.join(odir, 'barcodes.txt'), 'w') as f:
            f.write('\n'.join(adata.obs_names))
        np.savetxt(os.path.join(odir, 'spatial.tsv'), adata.obsm['spatial'], delimiter='\t')
        with open(os.path.join(odir, 'done.txt'), 'w') as f: f.write('ok')
        print(f"[{i+1}/{len(files)}] {tag} OK")
    except Exception as e:
        print(f"[{i+1}/{len(files)}] {tag} ERROR: {e}")

print("Done")
