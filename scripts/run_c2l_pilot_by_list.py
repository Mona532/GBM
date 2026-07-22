import argparse
from pathlib import Path
import warnings
import anndata as ad
import pandas as pd
import scanpy as sc
import cell2location
import torch

warnings.filterwarnings('ignore')
torch.set_float32_matmul_precision('medium')

p = argparse.ArgumentParser()
p.add_argument('--reference-dir', required=True)
p.add_argument('--sample-list', required=True)
p.add_argument('--outdir', required=True)
p.add_argument('--max-epochs', type=int, default=500)
args = p.parse_args()

reference_dir = Path(args.reference_dir)
sample_list = pd.read_csv(args.sample_list)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)
accel = 'gpu' if torch.cuda.is_available() else 'cpu'

sig = pd.read_csv(reference_dir / 'reference_signatures.csv', index_col=0)
genes = [x.strip() for x in (reference_dir / 'reference_genes.txt').read_text().splitlines() if x.strip()]

h5ad_dir = Path(r'E:/GBM/spatial_data_visium/spatial_data_visium/anndata')
vis_dir = Path(r'E:/GBM/ST_DATA/visium_all')
rows = []

for rec in sample_list.to_dict(orient='records'):
    sid = rec['sample']
    sample_out = outdir / sid
    sample_out.mkdir(parents=True, exist_ok=True)
    h5ad_path = h5ad_dir / f'{sid}.h5ad'
    vis_path = vis_dir / sid
    try:
        if h5ad_path.exists():
            adata = ad.read_h5ad(h5ad_path)
            adata.var_names_make_unique()
            if 'feature_types' in adata.var.columns:
                ge = adata[:, adata.var['feature_types'] == 'Gene Expression'].copy()
            else:
                ge = adata.copy()
        elif (vis_path / 'filtered_feature_bc_matrix.h5').exists():
            ge = sc.read_visium(vis_path)
            ge.var_names_make_unique()
            mt = ge.var_names.str.upper().str.startswith('MT-')
            ge = ge[:, ~mt].copy()
        else:
            raise FileNotFoundError(f'No input found for {sid}')
        shared = ge.var_names.intersection(sig.index).intersection(genes)
        ge = ge[:, shared].copy()
        sig_sub = sig.loc[shared]
        cell2location.models.Cell2location.setup_anndata(ge, batch_key=None)
        model = cell2location.models.Cell2location(ge, cell_state_df=sig_sub, N_cells_per_location=8, detection_alpha=200)
        model.train(max_epochs=args.max_epochs, batch_size=1024, train_size=1, accelerator=accel)
        ge = model.export_posterior(ge, sample_kwargs={'num_samples': 1000, 'batch_size': 1024})
        ct = list(ge.uns['mod']['factor_names'])
        for key, fname in [('q05_cell_abundance_w_sf', 'cell2loc_q05.csv'), ('means_cell_abundance_w_sf', 'cell2loc_mean.csv')]:
            arr = ge.obsm[key]
            vals = arr.values if hasattr(arr, 'values') else arr
            pd.DataFrame(vals, index=ge.obs_names, columns=ct).to_csv(sample_out / fname)
        ge.write(sample_out / 'cell2loc.h5ad', compression='gzip')
        rows.append({'sample': sid, 'status': 'done', 'spots': ge.n_obs, 'genes': len(shared)})
    except Exception as e:
        rows.append({'sample': sid, 'status': f'error: {str(e)[:120]}', 'spots': 0, 'genes': 0})

pd.DataFrame(rows).to_csv(outdir / 'pilot_summary.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
