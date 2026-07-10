"""Train 16-type reference (13 original + FDC + HEV + Tfh)."""
import anndata as ad, pandas as pd, numpy as np, cell2location, torch
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
torch.set_float32_matmul_precision("medium")

REF = Path(r"E:/GBM/results/c2l_consolidated_ref_tls16")
adata = ad.read_h5ad(REF / "reference_consolidated.h5ad")
print(f"Training: {adata.n_obs} cells, {adata.n_vars} genes, {adata.obs['c2l_label'].nunique()} types")

cell2location.models.RegressionModel.setup_anndata(adata, labels_key="c2l_label", batch_key="c2l_batch")
model = cell2location.models.RegressionModel(adata)
model.train(max_epochs=250, batch_size=2048, train_size=1, accelerator="gpu")

adata = model.export_posterior(adata, sample_kwargs={"num_samples":1000, "batch_size":2048})
means = adata.varm["means_per_cluster_mu_fg"]
fn = adata.uns["mod"]["factor_names"]
sig = pd.DataFrame(means.values, index=adata.var_names, columns=fn)
sig.to_csv(REF / "consolidated_signatures.csv")

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12)
forced = ["CR2","FCER2","CXCL13","FDCSP","ACKR1","SELE","CCL21","VCAM1","ICAM1","CHST4","MADCAM1",
          "BCL6","CXCR5","ICOS","PDCD1","IL21","SH2D1A","TOX2"]
selected = sorted(set(selected) | set(f for f in forced if f in adata.var_names))
(REF / "consolidated_genes.txt").write_text("\n".join(selected))
print(f"Done. Signatures: {sig.shape}, Selected genes: {len(selected)}")
