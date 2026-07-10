"""Assemble TLS reference v1 — backed mode, following working pattern."""
import anndata as ad, pandas as pd, numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix
import sys, warnings; warnings.filterwarnings("ignore")

OUT = Path(r"E:/GBM/results/reference_rebuild/tls_reference_v1")
OUT.mkdir(parents=True, exist_ok=True)

# ====== 1. Main GBM reference (backed mode, filter cells first) ======
print("Loading main GBM reference (backed)...")
main = ad.read_h5ad(r"E:/GBM/ST_DATA/GBM_space_snRNA/GBM_space_snRNA.h5ad", backed="r")

gbm_map = {
    "AC progenitor-like 1": "Glioma", "AC progenitor-like 2": "Glioma",
    "AC-gliosis-like 1": "Glioma", "NPC-neuronal-like 1": "Glioma",
    "NPC-neuronal-like 2": "Glioma", "NPC-neuronal-like 3": "Glioma",
    "NPC-neuronal-like 4": "Glioma", "Proliferative AC-OPC-like": "Glioma",
    "Proliferative nIPC-like": "Glioma", "Gliosis-like": "Glioma",
    "Hypoxic 1": "Glioma", "Hypoxic 2": "Glioma",
    "OPC-like 1": "Glial", "OPC-like 2": "Glial", "OPC-like 3": "Glial",
    "OPC-like 4": "Glial", "OPCs": "Glial", "OPC-NPC-like 1": "Glial",
    "OPC-NPC-like 2": "Glial", "OPC-NPC-like 3": "Glial",
    "OPC-neuronal-like": "Glial", "Oligodendrocytes 1": "Glial",
    "Oligodendrocytes 2": "Glial", "Ambiguous (oligo.)": "Glial",
    "Exc L2-3 IT": "Glial",
    "Endothelial (capillary)": "Vascular", "Endothelial (arteriole)": "Vascular",
    "Endothelial (venule)": "Vascular", "Endothelial (CNA-associated)": "Vascular",
    "Endothelial (Other)": "Vascular",
    "Pericytes 1": "Vascular", "Pericytes 2": "Vascular",
    "VLMC": "Vascular", "Ambiguous (vascular)": "Vascular",
    "Endothelial (capillary)": "Vascular", "Endothelial (arteriole)": "Vascular",
    "Endothelial (venule)": "Vascular", "Endothelial (CNA-associated)": "Vascular",
    "Endothelial (Other)": "Vascular",
    "Monocytes": "Macrophage",
    "T reg": "CD4_T", "Treg": "CD4_T", "CD4+ TEM cells": "CD4_T", "Naïve T cells": "CD4_T",
    "HSP-response T cells": "CD4_T", "IFN-response T cells": "CD4_T",
    "Proliferative T cells": "CD4_T",
    "CD8+ T cells": "CD8_T", "CD8+ T cells (cytotoxic)": "CD8_T",
    "NK cells 1": "NK", "NK cells 2": "NK",
    "B cells": "B", "Plasma cells": "Plasma",
    "Pro-inflammatory TAMs": "Macrophage", "Anti-inflammatory TAMs": "Macrophage",
    "Angiogenic TAMs": "Macrophage", "RTN1+ TAMs": "Macrophage",
    "Resident-TAMs": "Macrophage", "Resident BAM TAMs": "Macrophage",
    "Astrocyte-like TAMs": "Macrophage", "Interferon TAMs": "Macrophage",
    "Proliferative TAMs": "Macrophage", "Stress-response TAMs": "Macrophage",
    "Ambiguous (TAMs)": "Macrophage", "Dendritic cells": "cDC_or_mature_DC",
}

ct = main.obs["annotation_granular"].astype(str)
valid = ct.isin(gbm_map.keys())
main_idx = np.where(valid)[0]
print(f"Main filtered: {len(main_idx)} cells")

if main.raw:
    main_X = main.raw.X[main_idx, :]; main_var = main.raw.var_names
else:
    main_X = main.X[main_idx, :]; main_var = main.var_names

main_sub = ad.AnnData(X=main_X, obs=main.obs.iloc[main_idx].copy(),
                       var=pd.DataFrame(index=main_var))
main_sub.obs["ref_label"] = ct.iloc[main_idx].map(gbm_map).values
main_sub.obs["ref_batch"] = "gbm_main"
del main

# Cap
CAP = 3000; rng = np.random.RandomState(42)
keep = []
for lab in main_sub.obs["ref_label"].unique():
    idx = np.where(main_sub.obs["ref_label"] == lab)[0]
    if len(idx) > CAP: idx = rng.choice(idx, CAP, replace=False)
    keep.extend(idx)
main_sub = main_sub[sorted(keep)]
print(f"Main capped: {main_sub.n_obs} cells")
print(main_sub.obs["ref_label"].value_counts().to_string())

# ====== 2. Lymphoid organ reference ======
print("\nLoading lymphoid reference...")
lymph = ad.read_h5ad(r"E:/GBM/ST_DATA/lymphoid_organ_ref/sc.h5ad")
lymph_map = {
    "FDC": "FDC", "T_CD4+_TfH": "Tfh-like_CD4", "T_CD4+_TfH_GC": "Tfh-like_CD4",
    "T_Treg": "CD4_T", "B_GC_LZ": "B", "B_GC_DZ": "B", "B_GC_prePB": "B",
    "DC_cDC1": "cDC_or_mature_DC", "DC_cDC2": "cDC_or_mature_DC",
    "DC_CCR7+": "cDC_or_mature_DC", "DC_pDC": "cDC_or_mature_DC",
}
lymph.obs["ref_label"] = lymph.obs["Subset"].map(lymph_map)
lymph = lymph[lymph.obs["ref_label"].notna()].copy()
lymph.obs["ref_batch"] = "lymphoid_organ"
# Cap
keep = []
for lab in lymph.obs["ref_label"].unique():
    idx = np.where(lymph.obs["ref_label"] == lab)[0]
    if len(idx) > CAP: idx = rng.choice(idx, CAP, replace=False)
    keep.extend(idx)
lymph = lymph[sorted(keep)]
print(f"Lymphoid: {lymph.n_obs} cells")
print(lymph.obs["ref_label"].value_counts().to_string())

# ====== 3. HEV external ======
print("\nLoading HEV external...")
from scipy.io import mmread
hev_dir = Path(r"E:/GBM/results/reference_rebuild/tls_reference_v1/external_export")
hev_X = mmread(str(hev_dir / "stromal_counts.mtx")).tocsc()
hev_X.data = np.round(hev_X.data).astype(np.int32)
hev_labels = pd.read_csv(hev_dir / "stromal_labels.csv")
hev_genes = pd.read_csv(hev_dir / "stromal_genes.csv")["gene"].tolist()
hev_sub = ad.AnnData(X=hev_X.T.tocsr(), obs=hev_labels[["label"]],
                      var=pd.DataFrame(index=hev_genes))
hev_sub = hev_sub[hev_sub.obs["label"] == "HEV-like_endothelial"].copy()
hev_sub.obs["ref_label"] = "HEV-like_endothelial"
hev_sub.obs["ref_batch"] = "external_gbm"
print(f"HEV: {hev_sub.n_obs} cells, {hev_sub.n_vars} genes")

# ====== 4. ILC reference (backed, use raw) ======
print("\nLoading ILC reference (backed)...")
ilc = ad.read_h5ad(r"E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad", backed="r")
ilc_ct = ilc.obs["ilc_subtype"].astype(str)
ilc_keep = ilc_ct.isin(["ILC1","ILC2","ILC3"])
ilc_idx = np.where(ilc_keep)[0]
if ilc.raw:
    ilc_X = ilc.raw.X[ilc_idx, :]; ilc_var = ilc.raw.var_names
else:
    ilc_X = ilc.X[ilc_idx, :]; ilc_var = ilc.var_names
ilc_sub = ad.AnnData(X=ilc_X, obs=ilc.obs.iloc[ilc_idx].copy(),
                      var=pd.DataFrame(index=ilc_var))
ilc_sub.obs["ref_label"] = ilc_ct.iloc[ilc_idx].values
ilc_sub.obs["ref_batch"] = "ilc_ref"
del ilc
print(f"ILC: {ilc_sub.n_obs} cells")
print(ilc_sub.obs["ref_label"].value_counts().to_string())

# ====== 5. Merge on common genes ======
common = main_sub.var_names.intersection(lymph.var_names)
common = common.intersection(hev_sub.var_names).intersection(ilc_sub.var_names)
common = sorted(common)
print(f"\nCommon genes: {len(common)}")

merged = ad.concat([main_sub[:, common], lymph[:, common],
                     hev_sub[:, common], ilc_sub[:, common]], join="inner")

# ====== 6. Gene filter (sparse-safe, keep rare markers) ======
forced = pd.read_csv(r"E:/GBM/docs/tls_reference_forced_markers.csv")
forced_genes = forced["marker_gene"].dropna().unique()
# Sparse: count non-zero cells per gene
gene_n_cells = (merged.X > 0).sum(axis=0)
if hasattr(gene_n_cells, 'A1'): gene_n_cells = gene_n_cells.A1  # flatten
gene_ok = gene_n_cells >= max(10, merged.n_obs * 0.01)
# Force-keep forced markers
for g in forced_genes:
    if g in merged.var_names:
        gene_ok[list(merged.var_names).index(g)] = True
merged = merged[:, gene_ok]
print(f"After gene filter: {merged.n_vars} genes (forced={sum(gene_ok)})")

merged.write(OUT / "tls_reference_v1.h5ad", compression="gzip")
print(f"\nSaved: {merged.n_obs} cells, {merged.n_vars} genes, {merged.obs['ref_label'].nunique()} types")
print(merged.obs["ref_label"].value_counts().to_string())
