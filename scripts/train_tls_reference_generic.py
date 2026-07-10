import argparse
from pathlib import Path
import warnings
import anndata as ad
import cell2location
import pandas as pd
import torch

warnings.filterwarnings('ignore')
torch.set_float32_matmul_precision('medium')

p = argparse.ArgumentParser()
p.add_argument('--ref-h5ad', required=True)
p.add_argument('--outdir', required=True)
p.add_argument('--labels-key', default='ref_label')
p.add_argument('--batch-key', default='ref_source')
p.add_argument('--max-epochs', type=int, default=250)
args = p.parse_args()

ref_h5ad = Path(args.ref_h5ad)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(ref_h5ad)
print(f'Reference: {adata.n_obs} cells, {adata.n_vars} genes, {adata.obs[args.labels_key].nunique()} labels')
cell2location.models.RegressionModel.setup_anndata(adata, labels_key=args.labels_key, batch_key=args.batch_key)
model = cell2location.models.RegressionModel(adata)
model.train(max_epochs=args.max_epochs, batch_size=2048, train_size=1, accelerator='gpu' if torch.cuda.is_available() else 'cpu')
model.save(outdir / 'reference_model', overwrite=True)
adata = model.export_posterior(adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2048})
adata.write(outdir / 'reference_posterior.h5ad', compression='gzip')
means = adata.varm['means_per_cluster_mu_fg']
signatures = pd.DataFrame(means, index=adata.var_names, columns=adata.uns['mod']['factor_names'])
signatures.to_csv(outdir / 'reference_signatures.csv')
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12)
(outdir / 'reference_genes.txt').write_text('\n'.join(map(str, selected)))
print(f'Done. signatures={signatures.shape}, selected_genes={len(selected)}')
