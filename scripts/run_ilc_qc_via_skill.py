from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "ilc_qc_input_raw.h5ad"
OUTDIR = ROOT / "ilc_qc_only"
OUTPUT = OUTDIR / "ilc_qc_only_filtered_raw.h5ad"
SKILL_SCRIPT = Path(r"C:\Users\qingy\.agents\skills\h5ad-qc-preannotation\scripts\run_h5ad_qc_preannotation.py")


def load_skill_module():
    spec = importlib.util.spec_from_file_location("h5ad_qc_preannotation_skill", SKILL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load skill script: {SKILL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_args():
    class Args:
        create_min_features = 100
        nmads = 3.0
        hard_min_counts = None
        hard_min_features = 200
        hard_max_mt = 25.0
        hard_max_rb = 40.0
        hard_max_hb = 1.0
        high_feature_quantile = 0.99
        mt_prefix = "MT-"
        rb_regex = r"^RP[LS]"
        skip_scrublet = False
        fixed_doublet_rate = None
        doublet_rate_multiplier = 1.0
        doublet_min_cells = 50
        sim_doublet_ratio = 2.0
        stdev_doublet_rate = 0.02
        doublet_n_prin_comps = 30
        doublet_threshold = None
        doublet_calling_mode = "expected_adj_top_scores"
        random_seed = 42

    return Args()


def main() -> None:
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = True
    OUTDIR.mkdir(parents=True, exist_ok=True)

    skill = load_skill_module()
    skill.setup_logging()
    adata = sc.read_h5ad(INPUT)
    args = build_args()

    processed = []
    qc_stats = []
    doublet_stats = []
    for sample_id, idx in adata.obs.groupby("sample", observed=True).indices.items():
        sample_id = str(sample_id)
        skill.log_message(f"QC-only processing sample {sample_id}")
        sample = adata[list(idx), :].copy()
        sample.obs["sample"] = sample_id
        qc_sample, qc_df = skill.compute_qc_for_sample(sample, sample_id, args, "counts")
        qc_stats.append(qc_df)
        if qc_sample is None:
            skill.log_message(f"No cells left after QC for {sample_id}", "warning")
            continue
        singlets, dbl_df = skill.run_scrublet_for_sample(qc_sample, sample_id, args)
        processed.append(singlets)
        doublet_stats.append(dbl_df)

    if not processed:
        raise RuntimeError("No samples available after QC/doublet processing.")

    merged = ad.concat(processed, axis=0, join="outer", merge="same", index_unique=None)
    merged.obs_names_make_unique()
    merged.layers["counts"] = skill.get_counts_matrix(merged, "counts").copy()
    merged = skill.sanitize_for_h5ad_write(merged)
    merged.write_h5ad(OUTPUT, compression="gzip")

    pd.concat(qc_stats, ignore_index=True).to_csv(OUTDIR / "qc_statistics.tsv", sep="\t", index=False)
    pd.concat(doublet_stats, ignore_index=True).to_csv(OUTDIR / "doublet_statistics.tsv", sep="\t", index=False)

    print(f"Wrote {OUTPUT}")
    print(f"Cells={merged.n_obs} Genes={merged.n_vars}")


if __name__ == "__main__":
    main()
