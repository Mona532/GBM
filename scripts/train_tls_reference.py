"""Train cell2location reference model on TLS reference v1. Follows working pattern."""
import anndata as ad, pandas as pd, numpy as np, cell2location, torch
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

torch.set_float32_matmul_precision("medium")

REF = Path(r"E:/GBM/results/reference_rebuild/tls_reference_v1")
TRAINED = REF / "trained_reference"
TRAINED.mkdir(parents=True, exist_ok=True)

# Load reference
adata = ad.read_h5ad(REF / "tls_reference_v1.h5ad")
print(f"Reference: {adata.n_obs} cells, {adata.n_vars} genes, {adata.obs['ref_label'].nunique()} types")
print(adata.obs["ref_label"].value_counts().to_string())

# Setup + Train
cell2location.models.RegressionModel.setup_anndata(
    adata, labels_key="ref_label", batch_key="ref_batch"
)
model = cell2location.models.RegressionModel(adata)
model.train(max_epochs=250, batch_size=2500, train_size=1,
            accelerator="gpu" if torch.cuda.is_available() else "cpu")
model.save(str(TRAINED / "model"), overwrite=True)

# Export posterior
adata = model.export_posterior(adata, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
adata.write(TRAINED / "posterior.h5ad", compression="gzip")

# Export signatures
means = adata.varm["means_per_cluster_mu_fg"]
# Strip prefix from DataFrame column names
if hasattr(means, 'columns'):
    means.columns = [c.replace('means_per_cluster_mu_fg_', '') for c in means.columns]
signatures = means
signatures.to_csv(TRAINED / "signatures.csv")

# Gene selection — use low cutoffs to preserve rare cell type markers
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# Force-retain rare cell markers that got filtered out
forced_df = pd.read_csv(r"E:/GBM/docs/tls_reference_forced_markers.csv")
forced_genes = [g for g in forced_df["marker_gene"].unique() if g in adata.var_names]
lost = [g for g in forced_genes if g not in selected]
if lost:
    selected = list(selected) + lost
    print(f"Force-retained {len(lost)} rare markers: {lost}")

(TRAINED / "selected_genes.txt").write_text("\n".join(map(str, selected)))
print(f"Done. Signatures: {signatures.shape}, Selected genes: {len(selected)}")
