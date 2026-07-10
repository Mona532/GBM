from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


REGION_MAP = {
    "C": "Cortex",
    "T": "Tumor",
    "TC": "TumorCore",
    "TI": "TumorInfiltration",
}


def parse_sample_metadata(sample_name: str) -> dict[str, str]:
    clean = sample_name.lstrip("#")
    parts = clean.split("_")
    region_code = parts[-2]
    return {
        "sample_id": clean,
        "sample_dir_name": sample_name,
        "patient_id": parts[0],
        "region_code": region_code,
        "region_label": REGION_MAP.get(region_code, region_code),
        "idh_status": "IDHmutant" if "IDHMutant" in clean else "IDHunknown",
        "cohort": "UniversityClinicFreiburg",
    }


def load_one_sample(sample_dir: Path) -> tuple[ad.AnnData, dict[str, object]]:
    sample_name = sample_dir.name
    outs_dir = sample_dir / "outs"
    adata = sc.read_visium(outs_dir)
    meta = parse_sample_metadata(sample_name)

    adata.var_names_make_unique()
    adata.obs_names = [f"{meta['sample_id']}_{x}" for x in adata.obs_names]
    for key, value in meta.items():
        adata.obs[key] = value
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    adata.obs["patient_id"] = adata.obs["patient_id"].astype("category")
    adata.obs["region_code"] = adata.obs["region_code"].astype("category")
    adata.obs["region_label"] = adata.obs["region_label"].astype("category")
    adata.obs["idh_status"] = adata.obs["idh_status"].astype("category")
    adata.obs["cohort"] = adata.obs["cohort"].astype("category")

    # Give each sample a unique spatial key so concatenation can keep all images.
    old_key = next(iter(adata.uns["spatial"].keys()))
    adata.uns["spatial"] = {meta["sample_id"]: adata.uns["spatial"][old_key]}

    summary = {
        **meta,
        "n_spots": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "spatial_key": meta["sample_id"],
    }
    return adata, summary


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: python integrate_visium2_anndata.py <input_dir> <output_dir> [sample_regex]"
        )

    input_dir = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    sample_regex = sys.argv[3] if len(sys.argv) >= 4 else r"^#"
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(
        p for p in input_dir.iterdir() if p.is_dir() and re.search(sample_regex, p.name)
    )
    if not sample_dirs:
        raise SystemExit(f"No sample directories matched regex: {sample_regex}")

    adatas: list[ad.AnnData] = []
    summaries: list[dict[str, object]] = []

    for sample_dir in sample_dirs:
        print(f"[load] {sample_dir.name}")
        adata, summary = load_one_sample(sample_dir)
        adatas.append(adata)
        summaries.append(summary)

    merged = ad.concat(
        adatas,
        axis=0,
        join="outer",
        merge="same",
        uns_merge="unique",
        label="concat_sample",
        index_unique=None,
    )
    merged.obs["concat_sample"] = merged.obs["concat_sample"].astype("category")
    merged.uns["integration_info"] = {
        "source_dir": str(input_dir),
        "n_samples": len(adatas),
        "sample_ids": [x["sample_id"] for x in summaries],
    }

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(output_dir / "visium2_sample_summary.csv", index=False)
    summary_df.to_json(output_dir / "visium2_sample_summary.json", orient="records", indent=2)
    merged.write_h5ad(output_dir / "visium2_merged_raw.h5ad", compression="gzip")

    with open(output_dir / "integration_info.json", "w", encoding="utf-8") as handle:
        json.dump(merged.uns["integration_info"], handle, ensure_ascii=False, indent=2)

    print(f"[done] samples={len(adatas)} spots={merged.n_obs} genes={merged.n_vars}")


if __name__ == "__main__":
    main()
