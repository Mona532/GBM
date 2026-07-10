from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse

IN = Path(r'E:/GBM/results/reference_rebuild/shareable_unified_reference/shareable_unified_reference_input.h5ad')
OUT = Path(r'E:/GBM/results/reference_rebuild/shareable_unified_reference/audit')
OUT.mkdir(parents=True, exist_ok=True)

PANELS = {
    'B': ['MS4A1','CD79A','CD79B','CD74','HLA-DRA'],
    'Plasma': ['MZB1','JCHAIN','SDC1','XBP1','SSR4'],
    'CD8_T': ['CD3D','CD3E','TRAC','CD8A','NKG7'],
    'NK': ['NKG7','KLRD1','GNLY','PRF1','FCGR3A'],
    'helper_CD4': ['CD3D','CD3E','TRAC','IL7R','LTB'],
    'Tfh-like_CD4': ['ICOS','PDCD1','CXCR5','BCL6','SH2D1A'],
    'Treg': ['FOXP3','IL2RA','CTLA4','TIGIT','IKZF2'],
    'cDC_or_mature_DC': ['FCER1A','CD1C','HLA-DRA','CCR7','LAMP3'],
    'HEV-like_endothelial': ['ACKR1','SELE','VCAM1','SELP','PECAM1'],
    'vascular_endothelial': ['CLDN5','KDR','PECAM1','EMCN','VWF'],
    'Pericyte_VLMC': ['RGS5','CSPG4','PDGFRB','COL1A1','COL3A1'],
    'Macrophage': ['C1QA','C1QB','FCER1G','TYROBP','AIF1'],
    'Glioma': ['EGFR','SOX2','OLIG2','PDGFRA','NES'],
    'Glial': ['MBP','MOG','PLP1','AQP4','GFAP'],
    'ILC1': ['NKG7','KLRD1','IFNG','TBX21','TRAC'],
    'ILC2': ['IL7R','GATA3','KLRB1','LTB','IL1RL1'],
    'ILC3': ['IL7R','KIT','LTB','NCR2','TNFRSF18'],
}

adata = ad.read_h5ad(IN)
log_adata = adata.copy()
sc.pp.normalize_total(log_adata, target_sum=1e4)
sc.pp.log1p(log_adata)
X = log_adata.X.toarray() if sparse.issparse(log_adata.X) else np.asarray(log_adata.X)
df = pd.DataFrame(X, index=log_adata.obs_names, columns=log_adata.var_names)

rows = []
summary = []
counts = log_adata.obs['ref_label'].astype(str).value_counts().to_dict()
for label in sorted(log_adata.obs['ref_label'].astype(str).unique()):
    idx = log_adata.obs['ref_label'].astype(str).eq(label).values
    sub = df.loc[idx]
    genes = [g for g in PANELS.get(label, []) if g in sub.columns]
    mean_panel = float(sub[genes].mean().mean()) if genes else np.nan
    detect_panel = float((sub[genes] > 0).mean().mean()) if genes else np.nan
    warning = []
    if counts.get(label, 0) < 30:
        warning.append('very_low_n')
    elif counts.get(label, 0) < 80:
        warning.append('low_n')
    if genes and detect_panel < 0.15:
        warning.append('weak_panel_detection')
    if genes and mean_panel < 0.2:
        warning.append('weak_panel_expression')
    summary.append({'ref_label': label, 'n_cells': counts.get(label, 0), 'panel_genes_present': len(genes), 'panel_mean_lognorm': mean_panel, 'panel_detect_rate': detect_panel, 'warning': ';'.join(warning) if warning else ''})
    for gene in genes:
        rows.append({
            'ref_label': label,
            'gene': gene,
            'mean_lognorm': float(sub[gene].mean()),
            'detect_rate': float((sub[gene] > 0).mean()),
            'n_cells': counts.get(label, 0),
        })

marker_audit = pd.DataFrame(rows)
summary_df = pd.DataFrame(summary).sort_values(['warning','n_cells','ref_label'], ascending=[False, False, True])

marker_audit.to_csv(OUT / 'shareable_unified_reference_marker_audit.csv', index=False)
summary_df.to_csv(OUT / 'shareable_unified_reference_training_readiness.csv', index=False)
print(summary_df.to_string(index=False))
