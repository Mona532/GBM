"""Add FDC, Tfh, HEV to existing 13-type reference — TRUE union gene set."""
import anndata as ad, pandas as pd, numpy as np
from pathlib import Path
from scipy.io import mmread
from scipy.sparse import csr_matrix, lil_matrix
import warnings; warnings.filterwarnings("ignore")

OLD_REF = Path(r"E:/GBM/results/c2l_consolidated_ref")
OUT = Path(r"E:/GBM/results/c2l_consolidated_ref_tls16")
OUT.mkdir(parents=True, exist_ok=True)

# Load old reference
old = ad.read_h5ad(OLD_REF / "reference_consolidated.h5ad")
print(f"Old ref: {old.n_obs} cells, {old.n_vars} genes")

# Collect ALL genes from all sources
lymph = ad.read_h5ad(r"E:/GBM/ST_DATA/lymphoid_organ_ref/sc.h5ad")
ilc = ad.read_h5ad(r"E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad")
hev_genes = pd.read_csv(r"E:/GBM/results/reference_rebuild/tls_reference_v1/external_export/stromal_genes.csv")["gene"].tolist()
lymph_genes = list(lymph.raw.var_names) if lymph.raw else list(lymph.var_names)
ilc_genes = list(ilc.raw.var_names) if ilc.raw else list(ilc.var_names)

all_genes = sorted(set(old.var_names) | set(lymph_genes) | set(ilc_genes) | set(hev_genes))
print(f"Union genes: {len(all_genes)} (old={len(old.var_names)}, +lymph={len(lymph_genes)-len(set(old.var_names)&set(lymph_genes))}, +ilc={len(ilc_genes)-len(set(old.var_names)&set(ilc_genes))}, +hev={len(hev_genes)-len(set(old.var_names)&set(hev_genes))})")

gene_to_idx = {g: i for i, g in enumerate(all_genes)}

# ====== Helper: pad a source to the full gene set ======
def pad_to_union(src_X, src_genes, all_genes, dtype=np.float32):
    src_gi = {g: i for i, g in enumerate(src_genes)}
    X = lil_matrix((src_X.shape[0], len(all_genes)), dtype=dtype)
    for j, g in enumerate(all_genes):
        if g in src_gi:
            X[:, j] = src_X[:, src_gi[g]].toarray().ravel()
    return X.tocsr()

# ====== Old ref (already capped) ======
old_pad = pad_to_union(old.X, list(old.var_names), all_genes)
ref_labels = old.obs["c2l_label"].values.copy()
ref_batches = ["gbm_main"] * old.n_obs

# ====== FDC ======
fdc_cells = lymph[lymph.obs["Subset"] == "FDC"]
fdc_X = lymph.raw[fdc_cells.obs_names].X if lymph.raw else fdc_cells.X
fdc_pad = pad_to_union(fdc_X, list(lymph.raw.var_names if lymph.raw else fdc_cells.var_names), all_genes)
ref_labels = np.append(ref_labels, ["FDC"] * fdc_pad.shape[0])
ref_batches.extend(["lymphoid"] * fdc_pad.shape[0])
print(f"FDC: {fdc_pad.shape[0]} cells")

# ====== HEV ======
hev_mtx = mmread(str(Path(r"E:/GBM/results/reference_rebuild/tls_reference_v1/external_export") / "stromal_counts.mtx")).tocsc()
hev_labels = pd.read_csv(r"E:/GBM/results/reference_rebuild/tls_reference_v1/external_export/stromal_labels.csv")
hev_mask = (hev_labels["label"] == "HEV-like_endothelial").values
hev_X = hev_mtx[:, hev_mask].T.tocsr()
hev_pad = pad_to_union(hev_X, hev_genes, all_genes)
hev_pad.data = np.round(hev_pad.data).astype(np.int32)  # writeMM float fix
ref_labels = np.append(ref_labels, ["HEV-like_endothelial"] * hev_pad.shape[0])
ref_batches.extend(["external_gbm"] * hev_pad.shape[0])
print(f"HEV: {hev_pad.shape[0]} cells")

# ====== Tfh ======
tfh_cells = lymph[lymph.obs["Subset"].isin(["T_CD4+_TfH", "T_CD4+_TfH_GC"])]
tfh_X = lymph.raw[tfh_cells.obs_names].X if lymph.raw else tfh_cells.X
tfh_pad = pad_to_union(tfh_X, list(lymph.raw.var_names if lymph.raw else tfh_cells.var_names), all_genes)
ref_labels = np.append(ref_labels, ["Tfh-like_CD4"] * tfh_pad.shape[0])
ref_batches.extend(["lymphoid"] * tfh_pad.shape[0])
print(f"Tfh: {tfh_pad.shape[0]} cells")

# ====== Combine ======
X_all = csr_matrix(np.vstack([old_pad.toarray(), fdc_pad.toarray(), hev_pad.toarray(), tfh_pad.toarray()]))
combined = ad.AnnData(X=X_all, obs=pd.DataFrame({"c2l_label": ref_labels, "c2l_batch": ref_batches}),
                      var=pd.DataFrame(index=all_genes))
print(f"Combined: {combined.n_obs} cells, {combined.n_vars} genes")

# Cap at 3000
CAP = 3000; rng = np.random.RandomState(42)
keep = []
for lab in combined.obs["c2l_label"].unique():
    idx = np.where(combined.obs["c2l_label"] == lab)[0]
    if len(idx) > CAP: idx = rng.choice(idx, CAP, replace=False)
    keep.extend(idx)
combined = combined[sorted(keep)]
print(f"After capping: {combined.n_obs} cells")

# Verify integer
assert np.all(combined.X.data == combined.X.data.astype(int)), "Non-integer data!"
# Selected genes: all union genes (no filter during assembly, filter happens in training)
(OUT / "consolidated_genes.txt").write_text("\n".join(all_genes))
combined.write(OUT / "reference_consolidated.h5ad", compression="gzip")
print(f"Saved: {combined.n_obs} cells, {combined.n_vars} genes")
print(combined.obs["c2l_label"].value_counts().to_string())
