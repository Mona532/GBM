from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse, io

RNG = np.random.default_rng(42)
OUT = Path(r'E:/GBM/results/reference_rebuild/shareable_unified_reference')
OUT.mkdir(parents=True, exist_ok=True)

MAIN = Path(r'E:/GBM/ST_DATA/GBM_space_snRNA/GBM_space_snRNA.h5ad')
ILC = Path(r'E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad')
BLOCK_DIR = Path(r'E:/GBM/results/reference_rebuild/reference_blocks')
EXT_DIR = OUT / 'external_tls_support_export'

MAIN_LYMPH_MAP = {
    'B cells': 'B',
    'Plasma cells': 'Plasma',
    'CD8+ T cells': 'CD8_T',
    'CD8+ T cells (cytotoxic)': 'CD8_T',
    'NK cells 1': 'NK',
    'NK cells 2': 'NK',
}
GLIOMA_COARSE = {
    'AC-progenitor-like', 'AC-gliosis-like', 'Gliosis-like', 'Hypoxic',
    'NPC-neuronal-like', 'OPC-NPC-like', 'OPC-like', 'OPC-neuronal-like',
    'Proliferative', 'Undefined'
}
GLIAL_COARSE = {'Astrocytes', 'OPCs', 'Oligodendrocytes', 'Neurons (Exc)', 'Neurons (Inh)'}
MAIN_CAP = {'Glioma': 400, 'Glial': 300, 'Macrophage': 400, 'CD8_T': 250, 'NK': 180, 'B': 150, 'Plasma': 80}
HEV_MARKERS = ['ACKR1', 'SELE', 'CCL21', 'VCAM1', 'ICAM1', 'IL33', 'SELP']
TFH_MARKERS = ['BCL6', 'CXCR5', 'ICOS', 'PDCD1', 'IL21', 'SH2D1A', 'TOX2']
TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TIGIT']

def counts_layer(adata):
    if adata.raw is not None:
        return adata.raw.X, pd.Index(adata.raw.var_names)
    return adata.X, pd.Index(adata.var_names)

def to_lognorm_df(adata):
    b = adata.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    X = b.X.toarray() if sparse.issparse(b.X) else np.asarray(b.X)
    return pd.DataFrame(X, index=b.obs_names, columns=b.var_names)

def score_panel(df, genes):
    genes = [g for g in genes if g in df.columns]
    if not genes:
        return pd.Series(0.0, index=df.index), pd.Series(0, index=df.index)
    vals = df[genes]
    return vals.mean(axis=1), (vals > 0).sum(axis=1)

def finalize_part(adata, prefix):
    adata = adata.copy()
    adata.obs_names = [f'{prefix}::{x}' for x in adata.obs_names]
    return adata

parts = []
main = ad.read_h5ad(MAIN, backed='r')
main_X, main_var = counts_layer(main)
obs_main = main.obs[['annotation_coarse', 'annotation_granular', 'sample']].copy()
obs_main['orig_pos'] = np.arange(main.n_obs)
obs_main['ref_label'] = pd.Series(index=obs_main.index, dtype='object')
mask_lymph = obs_main['annotation_granular'].isin(MAIN_LYMPH_MAP)
obs_main.loc[mask_lymph, 'ref_label'] = obs_main.loc[mask_lymph, 'annotation_granular'].map(MAIN_LYMPH_MAP)
obs_main.loc[obs_main['annotation_coarse'].eq('Myeloid') & obs_main['ref_label'].isna(), 'ref_label'] = 'Macrophage'
obs_main.loc[obs_main['annotation_coarse'].isin(GLIOMA_COARSE) & obs_main['ref_label'].isna(), 'ref_label'] = 'Glioma'
obs_main.loc[obs_main['annotation_coarse'].isin(GLIAL_COARSE) & obs_main['ref_label'].isna(), 'ref_label'] = 'Glial'
obs_main = obs_main[obs_main['ref_label'].notna()].copy()
selected_orig = []
for lab, cap in MAIN_CAP.items():
    orig = obs_main.loc[obs_main['ref_label'].eq(lab), 'orig_pos'].to_numpy()
    if orig.size > cap:
        orig = RNG.choice(orig, size=cap, replace=False)
    selected_orig.extend(orig.tolist())
selected_orig = np.array(sorted(selected_orig), dtype=int)
main_sub_obs = obs_main.set_index('orig_pos').loc[selected_orig].copy()
main_sub = ad.AnnData(X=main_X[selected_orig, :], obs=main_sub_obs, var=pd.DataFrame(index=main_var))
main_sub.obs['ref_source'] = 'main_background_sampled'
parts.append(finalize_part(main_sub, 'main_bg'))

vasc = ad.read_h5ad(BLOCK_DIR / 'vascular_associated.h5ad')
vasc_df = to_lognorm_df(vasc)
hev_score, hev_detect = score_panel(vasc_df, HEV_MARKERS)
endo_mask = vasc.obs['annotation_granular'].astype(str).str.contains('Endothelial', regex=False)
ackr1_pos = (vasc_df['ACKR1'] > 0) if 'ACKR1' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
sele_pos = (vasc_df['SELE'] > 0) if 'SELE' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
selp_pos = (vasc_df['SELP'] > 0) if 'SELP' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
hev_cut = float(np.quantile(hev_score.loc[endo_mask], 0.92)) if endo_mask.sum() else np.inf
strict_hev = endo_mask & (hev_score >= hev_cut) & ((ackr1_pos | sele_pos | selp_pos).values)
vasc.obs['ref_label'] = np.where(endo_mask, 'vascular_endothelial', 'Pericyte_VLMC')
vasc.obs.loc[strict_hev, 'ref_label'] = 'HEV-like_endothelial'
vasc.obs['ref_source'] = 'main_vascular_block'
for lab, cap in [('HEV-like_endothelial', 99999), ('vascular_endothelial', 220), ('Pericyte_VLMC', 180)]:
    sub = vasc[vasc.obs['ref_label'].eq(lab)].copy()
    if sub.n_obs == 0:
        continue
    if sub.n_obs > cap:
        sub = sub[RNG.choice(np.arange(sub.n_obs), size=cap, replace=False)].copy()
    parts.append(finalize_part(sub, f'vasc_{lab}'))

cd4 = ad.read_h5ad(BLOCK_DIR / 'cd4_treg_context.h5ad')
cd4_df = to_lognorm_df(cd4)
tfh_score, tfh_detect = score_panel(cd4_df, TFH_MARKERS)
treg_score, treg_detect = score_panel(cd4_df, TREG_MARKERS)
seed_treg = cd4.obs['annotation_granular'].astype(str).eq('T reg')
treg_cut = float(np.quantile(treg_score, 0.90)) if cd4.n_obs else np.inf
cd4.obs['ref_label'] = 'helper_CD4'
cd4.obs.loc[seed_treg | ((treg_score >= treg_cut) & (treg_detect >= 2)), 'ref_label'] = 'Treg'
non_treg = cd4.obs['ref_label'].eq('helper_CD4')
tfh_cut = float(np.quantile(tfh_score.loc[non_treg], 0.90)) if non_treg.sum() else np.inf
cd4.obs.loc[non_treg & (tfh_score >= tfh_cut) & (tfh_detect >= 2), 'ref_label'] = 'Tfh-like_CD4'
cd4.obs['ref_source'] = 'main_cd4_block'
for lab, cap in [('Treg', 180), ('Tfh-like_CD4', 160), ('helper_CD4', 220)]:
    sub = cd4[cd4.obs['ref_label'].eq(lab)].copy()
    if sub.n_obs == 0:
        continue
    if sub.n_obs > cap:
        sub = sub[RNG.choice(np.arange(sub.n_obs), size=cap, replace=False)].copy()
    parts.append(finalize_part(sub, f'cd4_{lab}'))

cdc = ad.read_h5ad(BLOCK_DIR / 'dendritic_cells.h5ad')
cdc.obs['ref_label'] = 'cDC_or_mature_DC'
cdc.obs['ref_source'] = 'main_dendritic_block'
if cdc.n_obs > 220:
    cdc = cdc[RNG.choice(np.arange(cdc.n_obs), size=220, replace=False)].copy()
parts.append(finalize_part(cdc, 'cdc_main'))

ilc = ad.read_h5ad(ILC, backed='r')
ilc_mask = ilc.obs['ilc_subtype'].astype(str).isin(['ILC1', 'ILC2', 'ILC3']).values
ilc_idx = np.where(ilc_mask)[0]
ilc_X, ilc_var = counts_layer(ilc)
ilc_sub = ad.AnnData(X=ilc_X[ilc_idx, :], obs=ilc.obs.iloc[ilc_idx].copy(), var=pd.DataFrame(index=ilc_var))
ilc_sub.obs['ref_label'] = ilc_sub.obs['ilc_subtype'].astype(str)
ilc_sub.obs['ref_source'] = 'GBM_ilc_readonly'
parts.append(finalize_part(ilc_sub, 'ilc'))

ext_X = io.mmread(EXT_DIR / 'matrix.mtx').tocsr().T.tocsr()
ext_genes = pd.read_csv(EXT_DIR / 'features.tsv', sep='\t', header=None)[0].astype(str).tolist()
ext_cells = pd.read_csv(EXT_DIR / 'barcodes.tsv', sep='\t', header=None)[0].astype(str).tolist()
ext_meta = pd.read_csv(EXT_DIR / 'metadata.csv')
ext_meta.index = ext_cells
ext = ad.AnnData(X=ext_X, obs=ext_meta.loc[ext_cells].copy(), var=pd.DataFrame(index=ext_genes))
parts.append(finalize_part(ext, 'external'))

common = pd.Index(parts[0].var_names)
for p in parts[1:]:
    common = common.intersection(p.var_names)
merged = ad.concat([p[:, common].copy() for p in parts], join='inner', merge='same')
X = merged.X.toarray() if sparse.issparse(merged.X) else np.asarray(merged.X)
keep_gene = np.asarray((X > 0).sum(axis=0)).ravel() >= 5
merged = merged[:, keep_gene].copy()
merged.obs['ref_label'] = merged.obs['ref_label'].astype(str)
merged.obs['ref_source'] = merged.obs['ref_source'].astype(str)
label_counts = merged.obs['ref_label'].value_counts().rename_axis('ref_label').reset_index(name='n_cells')
source_counts = merged.obs.groupby(['ref_label', 'ref_source']).size().rename('n_cells').reset_index()
merged.write(OUT / 'shareable_unified_reference_input.h5ad', compression='gzip')
label_counts.to_csv(OUT / 'shareable_unified_reference_label_counts.csv', index=False)
source_counts.to_csv(OUT / 'shareable_unified_reference_source_counts.csv', index=False)
print('saved', OUT / 'shareable_unified_reference_input.h5ad')
print('shape', merged.shape)
print(label_counts.to_dict(orient='records'))
