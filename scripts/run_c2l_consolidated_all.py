"""Run cell2location with consolidated 13-type reference on all GBM samples"""
import pandas as pd, numpy as np, anndata as ad, cell2location, torch, scanpy as sc
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
torch.set_float32_matmul_precision("medium")

REF = Path(r"E:/GBM/results/c2l_consolidated_ref")
sig = pd.read_csv(REF / "consolidated_signatures.csv", index_col=0)
genes = [x.strip() for x in (REF / "consolidated_genes.txt").read_text().splitlines() if x.strip()]
accel = "gpu" if torch.cuda.is_available() else "cpu"
SKIP_DMG = {"GSE194329_DMG1","GSE194329_DMG2","GSE194329_DMG3","GSE194329_DMG4","GSE194329_DMG5"}
OUT = Path(r"E:/GBM/results/c2l_consolidated"); OUT.mkdir(parents=True, exist_ok=True)
SUMMARY = OUT / "summary.csv"

# === Dataset 1: dryad h5ad (97 samples) ===
h5ad_dir = Path(r"E:/GBM/spatial_data_visium/spatial_data_visium/anndata")
rows = []
for h5_path in sorted(h5ad_dir.glob("*.h5ad")):
    sid = h5_path.stem
    sample_out = OUT / sid; sample_out.mkdir(parents=True, exist_ok=True)
    if (sample_out / "cell2loc_q05.csv").exists():
        rows.append({"sample": sid, "spots": 0, "genes": 0, "status": "skipped"})
        continue
    try:
        print(f"[h5ad] {sid}")
        adata = ad.read_h5ad(h5_path); adata.var_names_make_unique()
        ge = adata[:, adata.var["feature_types"]=="Gene Expression"].copy()
        shared = ge.var_names.intersection(sig.index).intersection(genes)
        ge = ge[:, shared]; sig_sub = sig.loc[shared]
        cell2location.models.Cell2location.setup_anndata(ge, batch_key=None)
        model = cell2location.models.Cell2location(ge, cell_state_df=sig_sub, N_cells_per_location=8, detection_alpha=200)
        model.train(max_epochs=500, batch_size=1024, train_size=1, accelerator=accel)
        ge = model.export_posterior(ge, sample_kwargs={"num_samples":1000, "batch_size":1024})
        # Save CSV
        ct = list(ge.uns["mod"]["factor_names"])
        for key, fname in [("q05_cell_abundance_w_sf","cell2loc_q05.csv"), ("means_cell_abundance_w_sf","cell2loc_mean.csv")]:
            arr = ge.obsm[key]; vals = arr.values if hasattr(arr,"values") else arr
            pd.DataFrame(vals, index=ge.obs_names, columns=ct).to_csv(sample_out / fname)
        ge.write(sample_out / "cell2loc.h5ad", compression="gzip")
        rows.append({"sample": sid, "spots": ge.n_obs, "genes": len(shared), "status": "done"})
    except Exception as e:
        print(f"  FAIL: {e}")
        rows.append({"sample": sid, "spots": 0, "genes": 0, "status": f"error: {str(e)[:50]}"})

# === Dataset 2: visium_all (GBM tumor only, exclude DMG) ===
vis_dir = Path(r"E:/GBM/ST_DATA/visium_all")
for d in sorted(vis_dir.iterdir()):
    if not d.is_dir(): continue
    sid = d.name
    if sid in SKIP_DMG: continue
    if not (d / "filtered_feature_bc_matrix.h5").exists(): continue
    sample_out = OUT / sid; sample_out.mkdir(parents=True, exist_ok=True)
    if (sample_out / "cell2loc_q05.csv").exists():
        rows.append({"sample": sid, "spots": 0, "genes": 0, "status": "skipped"})
        continue
    try:
        print(f"[visium] {sid}")
        adata = sc.read_visium(d); adata.var_names_make_unique()
        mt = adata.var_names.str.upper().str.startswith("MT-")
        adata = adata[:, ~mt].copy()
        shared = adata.var_names.intersection(sig.index).intersection(genes)
        adata = adata[:, shared]; sig_sub = sig.loc[shared]
        cell2location.models.Cell2location.setup_anndata(adata, batch_key=None)
        model = cell2location.models.Cell2location(adata, cell_state_df=sig_sub, N_cells_per_location=8, detection_alpha=200)
        model.train(max_epochs=500, batch_size=1024, train_size=1, accelerator=accel)
        adata = model.export_posterior(adata, sample_kwargs={"num_samples":1000, "batch_size":1024})
        ct = list(adata.uns["mod"]["factor_names"])
        for key, fname in [("q05_cell_abundance_w_sf","cell2loc_q05.csv"), ("means_cell_abundance_w_sf","cell2loc_mean.csv")]:
            arr = adata.obsm[key]; vals = arr.values if hasattr(arr,"values") else arr
            pd.DataFrame(vals, index=adata.obs_names, columns=ct).to_csv(sample_out / fname)
        adata.write(sample_out / "cell2loc.h5ad", compression="gzip")
        rows.append({"sample": sid, "spots": adata.n_obs, "genes": len(shared), "status": "done"})
    except Exception as e:
        print(f"  FAIL: {e}")
        rows.append({"sample": sid, "spots": 0, "genes": 0, "status": f"error: {str(e)[:50]}"})

pd.DataFrame(rows).to_csv(SUMMARY, index=False)
print(f"Done. {sum(1 for r in rows if r['status']=='done')} done, {sum(1 for r in rows if r['status']=='skipped')} skipped")
