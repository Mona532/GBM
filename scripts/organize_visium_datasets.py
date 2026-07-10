"""Organize GSE237183 and GSE194329 Visium data into sc.read_visium() compatible directories"""
import os, shutil, tarfile, gzip
from pathlib import Path

ST_DATA = Path(r"E:\GBM\ST_DATA")
OUT = Path(r"E:\GBM\ST_DATA\visium_all")
OUT.mkdir(parents=True, exist_ok=True)

# ===== GSE237183: organize by sample =====
gse237183 = ST_DATA / "GSE237183_RAW"
samples_237183 = set()
for f in gse237183.glob("*_filtered_feature_bc_matrix.h5"):
    # GSM7596587_mgh258_filtered_feature_bc_matrix.h5 -> mgh258
    base = f.stem
    parts = base.split("_", 1)  # [GSM7596587, mgh258_filtered_feature_bc_matrix]
    sample = parts[1].replace("_filtered_feature_bc_matrix", "")
    samples_237183.add(sample)

print(f"GSE237183: {len(samples_237183)} samples")

for sample in sorted(samples_237183):
    out_dir = OUT / f"GSE237183_{sample}"
    out_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir = out_dir / "spatial"
    spatial_dir.mkdir(exist_ok=True)

    # Find all files for this sample
    for f in gse237183.iterdir():
        if sample in f.name:
            if "filtered_feature_bc_matrix.h5" in f.name:
                shutil.copy2(f, out_dir / "filtered_feature_bc_matrix.h5")
            elif "tissue_positions_list" in f.name:
                # gunzip and copy
                dest = spatial_dir / "tissue_positions_list.csv"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
            elif "scalefactors_json" in f.name:
                dest = spatial_dir / "scalefactors_json.json"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
            elif "tissue_lowres_image" in f.name:
                dest = spatial_dir / "tissue_lowres_image.png"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
            elif "tissue_hires_image" in f.name:
                dest = spatial_dir / "tissue_hires_image.png"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
            elif "detected_tissue_image" in f.name:
                dest = spatial_dir / "detected_tissue_image.jpg"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
            elif "aligned_fiducials" in f.name:
                dest = spatial_dir / "aligned_fiducials.jpg"
                with gzip.open(f, "rb") as gz:
                    with open(dest, "wb") as out_f:
                        out_f.write(gz.read())
    print(f"  {sample} -> {out_dir}")

# ===== GSE194329: extract tar.gz =====
gse194329 = ST_DATA / "GSE194329_RAW"
tars = sorted(gse194329.glob("*.tar.gz"))
print(f"\nGSE194329: {len(tars)} archives")

for tar_path in tars:
    # GSM5833533_GBM1_spaceranger_out.tar.gz -> GBM1
    name = tar_path.stem.replace(".tar", "").replace("_spaceranger_out", "")
    parts = name.split("_", 1)  # ["GSM5833533", "GBM1"]
    sample = parts[1] if len(parts) > 1 else name
    # Remove DMG prefix for consistent naming
    out_dir = OUT / f"GSE194329_{sample}"
    if out_dir.exists():
        print(f"  {sample} -> already exists, skip")
        continue

    print(f"  extracting {sample}...")
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=out_dir)
    # spaceranger output is inside 'outs/' subdirectory
    inner = out_dir / "outs"
    if inner.exists():
        # Move everything from outs/ up one level
        for item in inner.iterdir():
            shutil.move(str(item), str(out_dir / item.name))
        inner.rmdir()
    print(f"  {sample} -> {out_dir}")

print(f"\nDone! Total samples in {OUT}: {len(list(OUT.iterdir()))}")
