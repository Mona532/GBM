library(SeuratDisk)
library(Seurat)

DATA_DIR <- "E:/GBM/GSE237183_RAW"
OUT_DIR  <- "E:/GBM/results/cell2loc"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# First sample: mgh258
samp_gsm <- "GSM7596587"
samp_tag <- "mgh258"
file_prefix <- file.path(DATA_DIR, paste0(samp_gsm, "_", samp_tag))

# Setup Visium directory structure
tmp_dir <- file.path(OUT_DIR, "tmp", samp_gsm)
sp_dir <- file.path(tmp_dir, "spatial")
dir.create(sp_dir, recursive=TRUE, showWarnings=FALSE)

h5_file <- paste0(file_prefix, "_filtered_feature_bc_matrix.h5")
file.copy(h5_file, file.path(tmp_dir, "filtered_feature_bc_matrix.h5"), overwrite=TRUE)

# Unzip spatial support files
for(ext in c("tissue_positions_list.csv.gz","scalefactors_json.json.gz","tissue_lowres_image.png.gz")) {
  src <- paste0(file_prefix, "_", ext)
  destname <- gsub(".gz$", "", ext)
  dsn <- file.path(sp_dir, destname)
  if(file.exists(src)) R.utils::gunzip(src, destname=dsn, overwrite=TRUE, remove=FALSE)
}

# Load as Seurat Spatial
se <- Load10X_Spatial(data.dir=tmp_dir, assay="Spatial", filter.matrix=TRUE, slice=samp_tag)
cat(sprintf("Loaded %d spots x %d genes\n", ncol(se), nrow(se)))

# Save as h5ad
SaveH5Seurat(se, file.path(OUT_DIR, paste0(samp_tag, ".h5Seurat")), overwrite=TRUE)
Convert(file.path(OUT_DIR, paste0(samp_tag, ".h5Seurat")), dest="h5ad", overwrite=TRUE)
cat(sprintf("%s.h5ad saved\n", samp_tag))
unlink(tmp_dir, recursive=TRUE)
