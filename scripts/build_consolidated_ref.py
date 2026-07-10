"""Build consolidated reference: 10 backbone types + ILC1/2/3, per-subtype markers"""
import anndata as ad, pandas as pd, numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix

OUT = Path(r"E:/GBM/results/c2l_consolidated_ref"); OUT.mkdir(parents=True, exist_ok=True)

# === Main reference consolidation map ===
CONSOLIDATE = {
    "B cells": "B", "Plasma cells": "Plasma",
    "CD8+ T cells": "CD8_T", "CD8+ T cells (cytotoxic)": "CD8_T",
    "CD4+ T cells": "CD4_T", "Naïve T cells": "CD4_T", "HSP-response T cells": "CD4_T",
    "IFN-response T cells": "CD4_T", "Proliferative T cells": "CD4_T", "T reg": "CD4_T",
    "NK cells 1": "NK", "NK cells 2": "NK",
    "Dendritic cells": "Dendritic",
    "Monocytes": "Macrophage",
    "Angiogenic TAMs": "Macrophage", "Anti-inflammatory TAMs": "Macrophage",
    "Astrocyte-like TAMs": "Macrophage", "Interferon TAMs": "Macrophage",
    "Pro-inflammatory TAMs": "Macrophage", "Proliferative TAMs": "Macrophage",
    "RTN1+ TAMs": "Macrophage", "Resident BAM TAMs": "Macrophage",
    "Resident-TAMs": "Macrophage", "Stress-response TAMs": "Macrophage",
    "Ambiguous (TAMs)": "Macrophage",
    "AC progenitor-like": "Glioma", "AC-gliosis-like": "Glioma",
    "Gliosis-like": "Glioma", "NPC-neuronal-like": "Glioma",
    "Proliferative AC-OPC-like": "Glioma", "Proliferative NPC-OPC-like": "Glioma",
    "Proliferative nIPC-like": "Glioma",
    "Astrocytes": "Glial", "OPC-like": "Glial", "OPC-NPC-like": "Glial",
    "OPC-neuronal-like": "Glial", "OPCs": "Glial", "Oligodendrocytes": "Glial",
    "Neurons (Exc)": "Glial", "Neurons (Inh)": "Glial",
    "Endothelial": "Vascular", "Pericytes": "Vascular", "VLMC": "Vascular",
    "Vascular-associated": "Vascular",
}

# === Load main reference ===
# Read obs only first, filter, then load raw counts for filtered cells
main = ad.read_h5ad(r"E:/GBM/ST_DATA/GBM_space_snRNA/GBM_space_snRNA.h5ad", backed="r")
ct = main.obs["annotation_granular"].astype(str)
valid = ct.isin(CONSOLIDATE.keys())
filtered_idx = np.where(valid)[0]
print(f"Filtered: {len(filtered_idx)} / {main.n_obs} cells")

# Load raw counts only for filtered cells
if main.raw:
    raw_X = main.raw.X[filtered_idx, :]
    raw_var = main.raw.var_names
else:
    raw_X = main.X[filtered_idx, :]
    raw_var = main.var_names
main_sub = ad.AnnData(X=raw_X, obs=main.obs.iloc[filtered_idx].copy(), var=pd.DataFrame(index=raw_var))
main_sub.obs["c2l_label"] = ct[filtered_idx].map(CONSOLIDATE).values
main_sub.obs["c2l_batch"] = "main_ref"  # single batch for main reference
del main  # free memory
# Cap large labels at 1000 cells
rng = np.random.RandomState(42)
keep = []
for lab in main_sub.obs["c2l_label"].unique():
    idx = np.where(main_sub.obs["c2l_label"] == lab)[0]
    if len(idx) > 1000: idx = rng.choice(idx, 1000, replace=False)
    keep.extend(idx)
main_sub = main_sub[sorted(keep)]
print(f"Main: {main_sub.n_obs} cells, {main_sub.obs['c2l_label'].nunique()} labels")

# === ILC reference ===
ilc = ad.read_h5ad(r"E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad", backed="r")
ilc_ct = ilc.obs["ilc_subtype"].astype(str)
ilc_keep = ilc_ct.isin(["ILC1","ILC2","ILC3"])
ilc_idx = np.where(ilc_keep)[0]
print(f"ILC filtered: {len(ilc_idx)} / {ilc.n_obs} cells")
if ilc.raw:
    ilc_X = ilc.raw.X[ilc_idx, :]; ilc_var = ilc.raw.var_names
else:
    ilc_X = ilc.X[ilc_idx, :]; ilc_var = ilc.var_names
ilc_sub = ad.AnnData(X=ilc_X, obs=ilc.obs.iloc[ilc_idx].copy(), var=pd.DataFrame(index=ilc_var))
ilc_sub.obs["c2l_label"] = ilc_ct.iloc[ilc_idx].values
ilc_sub.obs["c2l_batch"] = "ilc_ref"  # single batch for ILC reference
del ilc

# === Merge ===
common = main_sub.var_names.intersection(ilc_sub.var_names)
merged = ad.concat([main_sub[:, common], ilc_sub[:, common]], join="inner")

# Per-subtype ILC markers — force-keep before gene filter
ILC_MARKERS = {
    "ILC1": ["TBX21","EOMES","CCL5","KLRF1","NCR1","CXCR3","XCL1","IFNG"],
    "ILC2": ["GATA3","BCL11B","IL1RL1","IL17RB","ICOS","HHEY1","IL9R","IL4"],
    "ILC3": ["RORC","BATF3","BATF","TOX2","IL1R1","IL23R","NRP1","LTA","IL22"],
}
all_markers = set(g for genes in ILC_MARKERS.values() for g in genes)
marker_in = [g for g in all_markers if g in merged.var_names]

# Gene filter — keep all marker genes
Xc = merged.X.toarray() if hasattr(merged.X, "toarray") else np.asarray(merged.X)
gene_ok = (Xc > 0).sum(axis=0) >= max(10, merged.n_obs * 0.01)
marker_idx = [list(merged.var_names).index(g) for g in marker_in]
gene_ok[marker_idx] = True  # force-keep markers
merged = merged[:, gene_ok]

merged.write(OUT / "reference_consolidated.h5ad", compression="gzip")
print(f"Saved: {merged.n_obs} cells, {merged.n_vars} genes, {merged.obs['c2l_label'].nunique()} types")
for lab in sorted(merged.obs["c2l_label"].unique()):
    print(f"  {lab}: {int((merged.obs['c2l_label']==lab).sum())}")
