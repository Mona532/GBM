# Convert GBM_5samples_filtered_clean.rds to h5ad
# Preserve raw counts + metadata
library(Seurat)
library(SeuratDisk)

# ── Input ──
rds_path <- "E:/GBM/GBM_DATA/5sample_integration/GBM_5samples_filtered_clean.rds"
out_dir  <- "E:/GBM/GBM_DATA/5sample_integration"
h5seurat_path <- file.path(out_dir, "GBM_5samples_filtered_clean.h5Seurat")
h5ad_path      <- file.path(out_dir, "GBM_5samples_filtered_clean.h5ad")

# ── Load ──
cat("Loading RDS...\n")
obj <- readRDS(rds_path)
cat(sprintf("Object: %d cells × %d features\n", ncol(obj), nrow(obj)))

# ── Inspect ──
cat("\n── Assays ──\n")
print(names(obj@assays))

cat("\n── Default assay ──\n")
print(DefaultAssay(obj))

cat("\n── Metadata columns ──\n")
print(head(colnames(obj@meta.data), 30))

cat("\n── Metadata dim ──\n")
print(dim(obj@meta.data))

# Check if counts slot exists in default assay
def_assay <- DefaultAssay(obj)
cat(sprintf("\n── Slots in %s assay ──\n", def_assay))
print(slotNames(obj@assays[[def_assay]]))

# Check if data (log-normalized) and counts (raw) are present
if ("counts" %in% slotNames(obj@assays[[def_assay]])) {
  cat(sprintf("\ncounts dim: %d × %d\n", nrow(obj@assays[[def_assay]]@counts), ncol(obj@assays[[def_assay]]@counts)))
}
if ("data" %in% slotNames(obj@assays[[def_assay]])) {
  cat(sprintf("data dim: %d × %d\n", nrow(obj@assays[[def_assay]]@data), ncol(obj@assays[[def_assay]]@data)))
}

# ── Ensure counts are in the default assay ──
# If counts are in a different assay (e.g., "RNA" has counts but default is "SCT" etc), switch
DefaultAssay(obj) <- def_assay

# ── Convert to h5ad ──
cat("\nSaving as h5Seurat...\n")
SaveH5Seurat(obj, h5seurat_path, overwrite = TRUE)

cat("Converting to h5ad...\n")
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

cat(sprintf("\nDone! Output: %s\n", h5ad_path))

# Clean up intermediate h5Seurat
file.remove(h5seurat_path)
cat("Intermediate h5Seurat cleaned up.\n")
