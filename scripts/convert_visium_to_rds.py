"""
Convert 142 GBM Visium h5ad files to Seurat-ready spaceranger directory layout.

Output per sample: results/visium_rds/{sample}/
  filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  spatial/
    tissue_hires_image.png
    tissue_lowres_image.png   (copy of hires)
    scalefactors_json.json
    tissue_positions.csv
"""
import h5py, numpy as np, pandas as pd, os, gzip, json, shutil
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from PIL import Image as PILImage

root   = "E:/GBM"
h5_dir = f"{root}/spatial_data_visium/spatial_data_visium/anndata_core"
out_dir = f"{root}/results/visium_rds"
os.makedirs(out_dir, exist_ok=True)

h5_files = sorted([f for f in os.listdir(h5_dir) if f.endswith('.h5ad')])
stats = []

for h5_name in h5_files:
    sid = h5_name.replace('.h5ad', '')
    sample_dir = f"{out_dir}/{sid}"
    fm_dir = f"{sample_dir}/filtered_feature_bc_matrix"
    sp_dir = f"{sample_dir}/spatial"
    os.makedirs(fm_dir, exist_ok=True)
    os.makedirs(sp_dir, exist_ok=True)

    try:
        with h5py.File(f"{h5_dir}/{h5_name}", 'r') as f:
            genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f['/var/_index'][:]]
            spots = [x.decode() if isinstance(x, bytes) else str(x) for x in f['/obs/_index'][:]]
            coords = f['/obsm/spatial'][:]
            n_genes, n_spots = len(genes), len(spots)

            # ---- filtered_feature_bc_matrix (gzipped MEX) ----
            d = f['/X/data'][:]; idx = f['/X/indices'][:]; iptr = f['/X/indptr'][:]
            X = csr_matrix((d, idx, iptr), shape=(n_spots, n_genes))
            mmwrite(f"{fm_dir}/matrix.mtx", X, field='integer')
            with open(f"{fm_dir}/matrix.mtx", 'rb') as fi, gzip.open(f"{fm_dir}/matrix.mtx.gz", 'wb') as fo:
                fo.writelines(fi)
            os.remove(f"{fm_dir}/matrix.mtx")
            with gzip.open(f"{fm_dir}/barcodes.tsv.gz", 'wt', encoding='utf-8') as fh:
                fh.write('\n'.join(spots))
            with gzip.open(f"{fm_dir}/features.tsv.gz", 'wt', encoding='utf-8') as fh:
                for g in genes:
                    fh.write(f"{g}\t{g}\tGene Expression\n")

            # ---- spatial files ----
            # hires image
            sp_key = list(f['/uns/spatial'].keys())[0]
            # hires image: handle both uint8 (AT series) and float32 0-1 (dryad/GSE)
            img_data = np.squeeze(f[f'/uns/spatial/{sp_key}/images/hires'][:])
            if img_data.dtype == np.float32 or img_data.dtype == np.float64:
                img_data = (img_data * 255).clip(0, 255).astype(np.uint8)
            else:
                img_data = np.clip(img_data, 0, 255).astype(np.uint8)
            if img_data.ndim == 2:
                PILImage.fromarray(img_data, 'L').save(f"{sp_dir}/tissue_hires_image.png")
            elif img_data.shape[2] == 4:
                PILImage.fromarray(img_data, 'RGBA').save(f"{sp_dir}/tissue_hires_image.png")
            else:
                PILImage.fromarray(img_data, 'RGB').save(f"{sp_dir}/tissue_hires_image.png")
            w, h = img_data.shape[:2]
            fs = os.path.getsize(f"{sp_dir}/tissue_hires_image.png") / 1e6
            stats.append({'sample': sid, 'width_px': w, 'height_px': h, 'file_size_mb': round(fs, 2)})

            # scalefactors (add lowres = hires for the copy)
            sf = {}
            for k, v in f[f'/uns/spatial/{sp_key}/scalefactors'].items():
                sf[k] = float(v[()]) if v.shape == () else v[()]
            if 'tissue_hires_scalef' in sf:
                sf.setdefault('tissue_lowres_scalef', sf['tissue_hires_scalef'])
            with open(f"{sp_dir}/scalefactors_json.json", 'w') as fh:
                json.dump(sf, fh)

            # tissue_positions.csv (Seurat-standard name)
            pd.DataFrame({
                'barcode': spots,
                'in_tissue': f['/obs/in_tissue'][:],
                'array_row': f['/obs/array_row'][:],
                'array_col': f['/obs/array_col'][:],
                'pxl_col_in_fullres': coords[:, 0],
                'pxl_row_in_fullres': coords[:, 1],
            }).to_csv(f"{sp_dir}/tissue_positions.csv", index=False)

        print(f"{sid}: {n_spots} spots, {w}x{h} hires")

    except Exception as e:
        print(f"ERROR {sid}: {e}")

# Summary table
if stats:
    df = pd.DataFrame(stats)
    df.to_csv(f"{out_dir}/image_summary.csv", index=False)
    print(f"\n{len(stats)}/{len(h5_files)} samples, "
          f"size {df.width_px.min()}x{df.height_px.min()}-{df.width_px.max()}x{df.height_px.max()}")
print("Done")
