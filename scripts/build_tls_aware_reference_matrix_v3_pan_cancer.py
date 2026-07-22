import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

RNG = np.random.default_rng(42)
OUT = Path(r'E:/GBM/results/reference_rebuild/tls_aware_reference_matrix_v3_pan_cancer')
OUT.mkdir(parents=True, exist_ok=True)

MAIN = Path(r'E:/GBM/ST_DATA/GBM_space_snRNA/GBM_space_snRNA.h5ad')
ILC = Path(r'E:/GBM/GBM_DATA/5sample_final/GBM_ilc.h5ad')
BLOCK_DIR = Path(r'E:/GBM/results/reference_rebuild/reference_blocks')

BG_MAP = {
    'B cells': 'B', 'Plasma cells': 'Plasma', 'CD8+ T cells': 'CD8_T', 'CD8+ T cells (cytotoxic)': 'CD8_T',
    'NK cells 1': 'NK', 'NK cells 2': 'NK', 'Monocytes': 'Macrophage', 'Angiogenic TAMs': 'Macrophage',
    'Anti-inflammatory TAMs': 'Macrophage', 'Astrocyte-like TAMs': 'Macrophage', 'Interferon TAMs': 'Macrophage',
    'Pro-inflammatory TAMs': 'Macrophage', 'Proliferative TAMs': 'Macrophage', 'RTN1+ TAMs': 'Macrophage',
    'Resident BAM TAMs': 'Macrophage', 'Resident-TAMs': 'Macrophage', 'Stress-response TAMs': 'Macrophage',
    'Ambiguous (TAMs)': 'Macrophage', 'AC progenitor-like': 'Glioma', 'AC-gliosis-like': 'Glioma',
    'Gliosis-like': 'Glioma', 'NPC-neuronal-like': 'Glioma', 'Proliferative AC-OPC-like': 'Glioma',
    'Proliferative NPC-OPC-like': 'Glioma', 'Proliferative nIPC-like': 'Glioma', 'Astrocytes': 'Glial',
    'OPC-like': 'Glial', 'OPC-NPC-like': 'Glial', 'OPC-neuronal-like': 'Glial', 'OPCs': 'Glial',
    'Oligodendrocytes': 'Glial', 'Neurons (Exc)': 'Glial', 'Neurons (Inh)': 'Glial'
}
BG_CAP = {'B':250,'Plasma':150,'CD8_T':400,'NK':300,'Macrophage':500,'Glioma':600,'Glial':500}
HEV_MARKERS = ['ACKR1', 'SELE', 'CCL21', 'VCAM1', 'ICAM1', 'IL33', 'SELP']
TFH_MARKERS = ['BCL6', 'CXCR5', 'ICOS', 'PDCD1', 'IL21', 'SH2D1A', 'TOX2']
TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TIGIT']
DC_MARKERS = ['FCER1A', 'CD1C', 'HLA-DRA', 'CCR7', 'LAMP3', 'CD86', 'CCL19']
PAN_FDC = ['CR2','FCER2','TNFSF13B','CXCL13','C7','CLU','VCAM1','ICAM1']
EXCL_ENDO = ['PECAM1','KDR','VWF','EMCN']
EXCL_PERI = ['RGS5','CSPG4','PDGFRB','MCAM','DES']


def lognorm_df(adata):
    b = adata.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    X = b.X.toarray() if hasattr(b.X, 'toarray') else np.asarray(b.X)
    return pd.DataFrame(X, index=b.obs_names, columns=b.var_names)


def score_panel(df, genes):
    gs = [g for g in genes if g in df.columns]
    if not gs:
        return pd.Series(0.0, index=df.index), pd.Series(0, index=df.index)
    vals = df[gs]
    return vals.mean(axis=1), (vals > 0).sum(axis=1)


def summarize_sources(obs_df):
    vc = obs_df['annotation_granular'].astype(str).value_counts()
    return '; '.join(f'{k}:{v}' for k, v in vc.head(8).items())

parts = []
qc_rows = []

main = ad.read_h5ad(MAIN, backed='r')
main_labels = main.obs['annotation_granular'].astype(str)
obs_main = main.obs[['annotation_granular', 'sample']].copy()
obs_main['target'] = main_labels.map(BG_MAP)
obs_main = obs_main[obs_main['target'].notna()].copy()
selected = []
for target, cap in BG_CAP.items():
    idx = np.where(obs_main['target'].values == target)[0]
    if len(idx) > cap:
        idx = RNG.choice(idx, size=cap, replace=False)
    selected.extend(idx.tolist())
selected = np.array(sorted(selected), dtype=int)
main_X = main.raw.X[selected, :] if main.raw is not None else main.X[selected, :]
main_var = main.raw.var_names if main.raw is not None else main.var_names
main_sub = ad.AnnData(X=main_X, obs=obs_main.iloc[selected].copy(), var=pd.DataFrame(index=main_var))
main_sub.obs['ref_label'] = main_sub.obs['target'].astype(str)
main_sub.obs['ref_source'] = 'main_background_sample'
parts.append(main_sub)
for lab, sdf in main_sub.obs.groupby('ref_label'):
    qc_rows.append({'ref_label': lab, 'n_cells': int(sdf.shape[0]), 'source': 'main_background_sample', 'source_labels': summarize_sources(sdf)})

vasc = ad.read_h5ad(BLOCK_DIR / 'vascular_associated.h5ad')
vasc_df = lognorm_df(vasc)
hev_score, hev_detect = score_panel(vasc_df, HEV_MARKERS)
fdc_score, fdc_det = score_panel(vasc_df, PAN_FDC)
endo_excl_score, endo_excl_det = score_panel(vasc_df, EXCL_ENDO)
peri_excl_score, peri_excl_det = score_panel(vasc_df, EXCL_PERI)
vasc.obs['hev_score'] = hev_score.values
vasc.obs['hev_detect'] = hev_detect.values
vasc.obs['fdc_pan_score'] = fdc_score.values
vasc.obs['fdc_pan_det'] = fdc_det.values
vasc.obs['endo_excl_det'] = endo_excl_det.values
vasc.obs['peri_excl_det'] = peri_excl_det.values
endo_mask = vasc.obs['annotation_granular'].astype(str).str.contains('Endothelial', regex=False)
ackr1_pos = vasc_df['ACKR1'] > 0 if 'ACKR1' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
sele_pos = vasc_df['SELE'] > 0 if 'SELE' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
selp_pos = vasc_df['SELP'] > 0 if 'SELP' in vasc_df.columns else pd.Series(False, index=vasc_df.index)
hev_cut = float(np.quantile(vasc.obs.loc[endo_mask, 'hev_score'], 0.92)) if endo_mask.sum() > 0 else np.inf
strict_hev = endo_mask & (vasc.obs['hev_score'] >= hev_cut) & ((ackr1_pos | sele_pos | selp_pos).values)
fdc_cut = float(np.quantile(vasc.obs['fdc_pan_score'], 0.995)) if vasc.n_obs > 0 else np.inf
strict_fdc = (
    (vasc.obs['fdc_pan_score'] >= fdc_cut) &
    (vasc.obs['fdc_pan_det'] >= 2) &
    (vasc.obs['endo_excl_det'] <= 1) &
    (vasc.obs['peri_excl_det'] <= 2) &
    (~endo_mask)
)
vasc.obs['ref_label'] = 'perivascular_stromal'
vasc.obs.loc[endo_mask, 'ref_label'] = 'vascular_endothelial'
vasc.obs.loc[strict_hev, 'ref_label'] = 'HEV-like_endothelial'
vasc.obs.loc[strict_fdc, 'ref_label'] = 'FDC_like_pan_cancer_guided'
vasc.obs['ref_source'] = 'vascular_block_plus_pan_cancer_FDC'
for lab in ['FDC_like_pan_cancer_guided', 'HEV-like_endothelial', 'vascular_endothelial', 'perivascular_stromal']:
    sub = vasc[vasc.obs['ref_label'] == lab].copy()
    if sub.n_obs == 0:
        continue
    if lab in {'vascular_endothelial', 'perivascular_stromal'} and sub.n_obs > 350:
        sub = sub[RNG.choice(np.arange(sub.n_obs), size=350, replace=False)].copy()
    parts.append(sub)
    qc_rows.append({'ref_label': lab, 'n_cells': int(sub.n_obs), 'source': 'vascular_block_plus_pan_cancer_FDC', 'source_labels': summarize_sources(sub.obs), 'mean_hev_score': float(sub.obs['hev_score'].mean()) if 'hev_score' in sub.obs else np.nan, 'mean_fdc_pan_score': float(sub.obs['fdc_pan_score'].mean()) if 'fdc_pan_score' in sub.obs else np.nan})

cd4 = ad.read_h5ad(BLOCK_DIR / 'cd4_treg_context.h5ad')
cd4_df = lognorm_df(cd4)
tfh_score, tfh_detect = score_panel(cd4_df, TFH_MARKERS)
treg_score, treg_detect = score_panel(cd4_df, TREG_MARKERS)
cd4.obs['tfh_score'] = tfh_score.values
cd4.obs['tfh_detect'] = tfh_detect.values
cd4.obs['treg_score'] = treg_score.values
cd4.obs['treg_detect'] = treg_detect.values
cd4.obs['ref_label'] = 'CD4_T'
seed_treg = cd4.obs['annotation_granular'].astype(str).eq('T reg')
treg_cut = float(np.quantile(cd4.obs['treg_score'], 0.90))
cd4.obs.loc[seed_treg | ((cd4.obs['treg_score'] >= treg_cut) & (cd4.obs['treg_detect'] >= 2)), 'ref_label'] = 'Treg'
non_treg = cd4.obs['ref_label'].eq('CD4_T')
tfh_cut = float(np.quantile(cd4.obs.loc[non_treg, 'tfh_score'], 0.90)) if non_treg.sum() > 0 else np.inf
cd4.obs.loc[non_treg & (cd4.obs['tfh_score'] >= tfh_cut) & (cd4.obs['tfh_detect'] >= 2), 'ref_label'] = 'Tfh-like_CD4'
cd4.obs['ref_source'] = 'cd4_treg_block'
for lab in ['Treg', 'Tfh-like_CD4', 'CD4_T']:
    sub = cd4[cd4.obs['ref_label'] == lab].copy()
    if sub.n_obs == 0:
        continue
    if lab == 'CD4_T' and sub.n_obs > 350:
        sub = sub[RNG.choice(np.arange(sub.n_obs), size=350, replace=False)].copy()
    parts.append(sub)
    qc_rows.append({'ref_label': lab, 'n_cells': int(sub.n_obs), 'source': 'cd4_treg_block', 'source_labels': summarize_sources(sub.obs), 'mean_tfh_score': float(sub.obs['tfh_score'].mean()), 'mean_treg_score': float(sub.obs['treg_score'].mean())})

cdc = ad.read_h5ad(BLOCK_DIR / 'dendritic_cells.h5ad')
cdc_df = lognorm_df(cdc)
dc_score, dc_detect = score_panel(cdc_df, DC_MARKERS)
cdc.obs['dc_score'] = dc_score.values
cdc.obs['dc_detect'] = dc_detect.values
cdc.obs['ref_label'] = 'cDC_or_mature_DC'
cdc.obs['ref_source'] = 'dendritic_block'
parts.append(cdc)
qc_rows.append({'ref_label': 'cDC_or_mature_DC', 'n_cells': int(cdc.n_obs), 'source': 'dendritic_block', 'source_labels': summarize_sources(cdc.obs), 'mean_dc_score': float(cdc.obs['dc_score'].mean())})

ilc = ad.read_h5ad(ILC, backed='r')
ilc_idx = np.where(ilc.obs['ilc_subtype'].astype(str).isin(['ILC1', 'ILC2', 'ILC3']).values)[0]
ilc_X = ilc.raw.X[ilc_idx, :] if ilc.raw is not None else ilc.X[ilc_idx, :]
ilc_var = ilc.raw.var_names if ilc.raw is not None else ilc.var_names
ilc_sub = ad.AnnData(X=ilc_X, obs=ilc.obs.iloc[ilc_idx].copy(), var=pd.DataFrame(index=ilc_var))
ilc_sub.obs['ref_label'] = ilc_sub.obs['ilc_subtype'].astype(str)
ilc_sub.obs['ref_source'] = 'GBM_ilc_readonly'
parts.append(ilc_sub)
for lab, sdf in ilc_sub.obs.groupby('ref_label'):
    qc_rows.append({'ref_label': lab, 'n_cells': int(sdf.shape[0]), 'source': 'GBM_ilc_readonly', 'source_labels': lab})

common = pd.Index(parts[0].var_names)
for p in parts[1:]:
    common = common.intersection(p.var_names)
merged = ad.concat([p[:, common].copy() for p in parts], join='inner', merge='same')
X = merged.X.toarray() if hasattr(merged.X, 'toarray') else np.asarray(merged.X)
keep_gene = np.asarray((X > 0).sum(axis=0)).ravel() >= 5
merged = merged[:, keep_gene].copy()
X = merged.X.toarray() if hasattr(merged.X, 'toarray') else np.asarray(merged.X)
label_order = sorted(merged.obs['ref_label'].astype(str).unique())
gene_means = []
for lab in label_order:
    idx = np.where(merged.obs['ref_label'].astype(str).values == lab)[0]
    gene_means.append(X[idx].mean(axis=0))
gene_means = pd.DataFrame(np.vstack(gene_means).T, index=merged.var_names, columns=label_order)
label_counts = merged.obs['ref_label'].astype(str).value_counts().rename_axis('ref_label').reset_index(name='n_cells')
qc = pd.DataFrame(qc_rows).drop_duplicates(subset=['ref_label', 'source', 'n_cells'])
qc = qc.merge(label_counts, on='ref_label', how='left', suffixes=('_premerge', '_final'))

merged.write(OUT / 'tls_aware_reference_cells_v3_pan_cancer.h5ad', compression='gzip')
gene_means.to_csv(OUT / 'tls_aware_reference_gene_means_v3_pan_cancer.csv.gz')
label_counts.to_csv(OUT / 'tls_aware_reference_label_counts_v3_pan_cancer.csv', index=False)
qc.to_csv(OUT / 'tls_aware_reference_qc_summary_v3_pan_cancer.csv', index=False)

print('saved_h5ad', OUT / 'tls_aware_reference_cells_v3_pan_cancer.h5ad')
print('shape', merged.shape)
print('labels', label_counts.to_dict(orient='records'))
