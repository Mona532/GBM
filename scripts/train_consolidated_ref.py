"""Train cell2location reference model on consolidated reference (13 types)"""
import anndata as ad, pandas as pd, numpy as np, cell2location, torch
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

torch.set_float32_matmul_precision("medium")
REF = Path(r"E:/GBM/results/c2l_consolidated_ref")
OUT = REF

# Load reference
adata = ad.read_h5ad(REF / "reference_consolidated.h5ad")
print(f"Reference: {adata.n_obs} cells, {adata.n_vars} genes, {adata.obs['c2l_label'].nunique()} types")

# Train RegressionModel
cell2location.models.RegressionModel.setup_anndata(
    adata, labels_key="c2l_label", batch_key="c2l_batch"
)
model = cell2location.models.RegressionModel(adata)
model.train(max_epochs=250, batch_size=2048, train_size=1, accelerator="gpu" if torch.cuda.is_available() else "cpu")
model.save(OUT / "consolidated_model", overwrite=True)

# Export posterior
adata = model.export_posterior(adata, sample_kwargs={"num_samples": 1000, "batch_size": 2048})
adata.write(OUT / "consolidated_posterior.h5ad", compression="gzip")

# Export signatures
means = adata.varm["means_per_cluster_mu_fg"]
signatures = pd.DataFrame(
    means, index=adata.var_names,
    columns=adata.uns["mod"]["factor_names"]
)
signatures.to_csv(OUT / "consolidated_signatures.csv")

# Gene selection
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12)
(OUT / "consolidated_genes.txt").write_text("\n".join(map(str, selected)))
print(f"Done. Signatures: {signatures.shape}, Selected genes: {len(selected)}")
