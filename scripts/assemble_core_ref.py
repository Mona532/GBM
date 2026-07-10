"""Assemble 17-type reference. Gene dedup = same name => sum counts."""
import scanpy as sc, pandas as pd, numpy as np
from pathlib import Path
from scipy.io import mmread
from scipy.sparse import csr_matrix, lil_matrix, hstack, vstack
import warnings; warnings.filterwarnings("ignore")

OUT = Path(r"E:/GBM/results/reference_rebuild/core_ref_v1"); OUT.mkdir(parents=True, exist_ok=True)
CAP = 1000; rng = np.random.RandomState(42)

def dedup_genes(adata):
    """Collapse duplicate gene names by summing counts (QC skill pattern)."""
    gmap = {}
    for i, g in enumerate(adata.var_names):
        gmap.setdefault(g, []).append(i)
    if len(gmap) == adata.n_vars:
        return adata  # no duplicates
    X_new = hstack([adata.X[:, idxs].sum(axis=1) for idxs in gmap.values()])
    return sc.AnnData(X=X_new, obs=adata.obs, var=pd.DataFrame(index=list(gmap.keys())))

# ====== 1. CORE_GBM ======
print("CORE...")
core = sc.read_h5ad(r"E:/GBM/REF/CORE_GBM.h5ad", backed="r")
core.var_names = core.var["feature_name"].astype(str).values
valid = (core.var_names!="nan") & ~pd.Series(core.var_names).str.startswith("ENSG").values
core = core[:, valid].to_memory()
rawX = core.raw[:, valid].X.copy()
rawVar = core.raw[:, valid].var.copy()
core.X = csr_matrix(rawX); core.var = rawVar; del core.raw, rawX, rawVar
core = dedup_genes(core)

label_map = {"malignant cell":"Glioma","macrophage":"Macrophage","microglial cell":"Glial",
    "oligodendrocyte":"Glial","oligodendrocyte precursor cell":"Glial","radial glial cell":"Glial",
    "astrocyte":"Glial","neuron":"Glial","monocyte":"Monocyte","dendritic cell":"Dendritic",
    "natural killer cell":"NK","mural cell":"Vascular","endothelial cell":"Vascular",
    "B cell":"B","plasma cell":"Plasma","mast cell":"Macrophage"}
labels = core.obs["cell_type"].map(label_map).astype(str)
t = core.obs["cell_type"]=="mature T cell"
t_map = {"CD8 cytotoxic":"CD8_T","CD8 EM":"CD8_T","CD8 NK sig":"CD8_T",
    "CD4 INF":"CD4_T","CD4 rest":"CD4_T","Reg T":"CD4_T","Stress sig":"CD4_T","Prolif T":"CD4_T"}
labels[t] = core.obs.loc[t,"annotation_level_4"].map(t_map).fillna("CD4_T")
keep = []
for lab in labels.unique():
    idx = np.where(labels==lab)[0]; rng.shuffle(idx)
    if len(idx)>CAP: idx = idx[:CAP]
    keep.extend(idx)
core = core[sorted(keep)].copy()
core.obs = pd.DataFrame({"ref_label":labels.iloc[sorted(keep)].values, "ref_batch":"core_gbm"})
print(f"CORE: {core.n_obs}c, {core.n_vars}g")

# ====== 2. ILC ======
ilc = sc.read_h5ad(r"E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad")
ilc.obs["ref_label"] = ilc.obs["ilc_subtype"].map({"ILC1":"ILC1","ILC2":"ILC2","ILC3":"ILC3"})
ilc = ilc[ilc.obs["ref_label"].notna()].copy()
ilc.X = csr_matrix(ilc.raw.X.copy()); ilc.var = ilc.raw.var.copy(); del ilc.raw
ilc = dedup_genes(ilc)
ilc.obs = pd.DataFrame({"ref_label":ilc.obs["ref_label"].values, "ref_batch":"ilc_ref"})
print(f"ILC: {ilc.n_obs}")

# ====== 3. Lymphoid ======
lymph = sc.read_h5ad(r"E:/GBM/ST_DATA/lymphoid_organ_ref/sc.h5ad")
lymph.X = csr_matrix(lymph.raw.X.copy()); lymph.var = lymph.raw.var.copy(); del lymph.raw
lymph = dedup_genes(lymph)
fdc = lymph[lymph.obs["Subset"]=="FDC"].copy(); fdc.obs = pd.DataFrame({"ref_label":"FDC","ref_batch":"lymphoid"}, index=fdc.obs_names)
tfh = lymph[lymph.obs["Subset"].isin(["T_CD4+_TfH","T_CD4+_TfH_GC"])].copy()
tfh = tfh[rng.choice(tfh.n_obs, min(CAP,tfh.n_obs), replace=False)].copy()
tfh.obs = pd.DataFrame({"ref_label":"Tfh-like_CD4","ref_batch":"lymphoid"}, index=tfh.obs_names)
del lymph; print(f"FDC: {fdc.n_obs}, Tfh: {tfh.n_obs}")

# ====== 4. HEV ======
hev_dir = Path(r"E:/GBM/results/reference_rebuild/tls_reference_v1/external_export")
hev_X = mmread(str(hev_dir/"stromal_counts.mtx")).tocsc()
hev_lab = pd.read_csv(hev_dir/"stromal_labels.csv")
hev_gn = np.array(pd.read_csv(hev_dir/"stromal_genes.csv")["gene"].tolist())
hev_m = hev_lab["label"]=="HEV-like_endothelial"
hev = sc.AnnData(X=hev_X[:, hev_m].T.tocsr(), obs=pd.DataFrame({"ref_label":"HEV-like_endothelial","ref_batch":"external_gbm"}, index=range(hev_m.sum())), var=pd.DataFrame(index=hev_gn))
hev.X.data = np.round(hev.X.data).astype(np.int32)
hev = dedup_genes(hev); print(f"HEV: {hev.n_obs}")

# ====== 5. Concat outer ======
print("Concatenating...")
combined = sc.concat({"core":core,"ilc":ilc,"fdc":fdc,"tfh":tfh,"hev":hev}, join="outer", label="batch", index_unique="-")
all_var = pd.concat([v.var for v in [core,ilc,fdc,tfh,hev]], join="outer")
all_var = all_var[~all_var.index.duplicated()]
combined.var = all_var.loc[combined.var_names]
combined.X.data = np.nan_to_num(combined.X.data, nan=0).astype(np.int32)
combined.var_names = combined.var.index.values  # reset from var

# Final global dedup (same gene from diff sources)
combined = dedup_genes(combined)
suff = sum(1 for g in combined.var_names if g.endswith("-1"))
combined.write(OUT/"reference_consolidated.h5ad", compression="gzip")
print(f"Saved: {combined.n_obs}c, {combined.n_vars}g, {combined.obs['ref_label'].nunique()} types, suffix={suff}")
print(combined.obs["ref_label"].value_counts().to_string())
