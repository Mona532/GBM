"""Train 17-type reference + test one sample."""
import anndata as ad, pandas as pd, numpy as np, cell2location, torch, scanpy as sc, os, warnings
warnings.filterwarnings("ignore"); torch.set_float32_matmul_precision("medium")

REF = "E:/GBM/results/reference_rebuild/core_ref_v1"
os.makedirs(f"{REF}/trained", exist_ok=True)

# === Train ===
print("Training reference...")
adata = ad.read_h5ad(f"{REF}/reference_consolidated.h5ad")
adata.obs["c2l_label"] = adata.obs["ref_label"]
adata.obs["c2l_batch"] = adata.obs["ref_batch"]
cell2location.models.RegressionModel.setup_anndata(adata, labels_key="c2l_label", batch_key="c2l_batch")
model = cell2location.models.RegressionModel(adata)
model.train(max_epochs=250, batch_size=2048, train_size=1, accelerator="gpu")
adata = model.export_posterior(adata, sample_kwargs={"num_samples":1000, "batch_size":2048})
means = adata.varm["means_per_cluster_mu_fg"]
sig = pd.DataFrame(means.values, index=adata.var_names, columns=adata.uns["mod"]["factor_names"])
sig.to_csv(f"{REF}/signatures.csv")
from cell2location.utils.filtering import filter_genes
sel = list(filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12))
forced = ["SOX2","EGFR","OLIG2","MBP","PLP1","PECAM1","CDH5","CD68","CD4","CD3D",
          "CXCR5","PDCD1","BCL6","CD8A","CD8B","NKG7","PRF1","KLRD1",
          "MS4A1","CD79A","MZB1","JCHAIN","FCER1A","TBX21","GATA3","RORC",
          "CR2","FDCSP","ACKR1","SELE","CCL21"]
sel = sorted(set(sel) | {g for g in forced if g in adata.var_names})
open(f"{REF}/selected_genes.txt", "w").write("\n".join(sel))
print(f"Trained: {sig.shape[0]}g x {sig.shape[1]} types, selected {len(sel)} genes")

# === Test one sample ===
print("\nTesting one sample...")
sig = pd.read_csv(f"{REF}/signatures.csv", index_col=0)
genes = open(f"{REF}/selected_genes.txt").read().splitlines()
a = sc.read_h5ad("E:/GBM/spatial_data_visium/spatial_data_visium/anndata/AT10-BRA-5-FO-1_1.h5ad")
a.var_names_make_unique()
ge = a[:, a.var["feature_types"]=="Gene Expression"].copy()
shared = ge.var_names.intersection(sig.index).intersection(genes)
ge = ge[:, shared]; sig_sub = sig.loc[shared]
print(f"Shared genes: {len(shared)}")
cell2location.models.Cell2location.setup_anndata(ge, batch_key=None)
model = cell2location.models.Cell2location(ge, cell_state_df=sig_sub, N_cells_per_location=8, detection_alpha=200)
model.train(max_epochs=500, batch_size=2048, train_size=1, accelerator="gpu")
ge = model.export_posterior(ge, sample_kwargs={"num_samples":1000, "batch_size":2048})
q05 = ge.obsm["q05_cell_abundance_w_sf"]; vals = q05.values; ct = ge.uns["mod"]["factor_names"]
for i, c in enumerate(ct):
    print(f"{c:25s}: {100*(vals[:,i]>0).mean():.0f}% mean={vals[:,i].mean():.3f}")
print("Done")
