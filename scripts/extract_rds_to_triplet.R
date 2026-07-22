# Extract GBM_5samples_filtered_clean.rds → triplet files (counts + metadata + features)
library(Seurat)
library(Matrix)

rds_path <- "E:/GBM/GBM_DATA/5sample_integration/GBM_5samples_filtered_clean.rds"
out_dir  <- "E:/GBM/GBM_DATA/5sample_integration/triplet"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading RDS...\n")
obj <- readRDS(rds_path)
cat(sprintf("Object: %d cells × %d features\n", ncol(obj), nrow(obj)))

DefaultAssay(obj) <- "RNA"

# ── 1. Extract raw counts matrix (sparse) ──
cat("Extracting raw counts matrix...\n")
counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
cat(sprintf("Counts dim: %d × %d, sparsity: %.2f%%\n",
            nrow(counts), ncol(counts),
            100 * (1 - length(counts@x) / (nrow(counts) * ncol(counts)))))

# Save as MatrixMarket
writeMM(counts, file.path(out_dir, "counts.mtx"))
cat("  → counts.mtx\n")

# ── 2. Barcodes ──
write.table(colnames(counts), file.path(out_dir, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("  → barcodes.tsv\n")

# ── 3. Features (gene names) ──
feat_df <- data.frame(
  gene_id = rownames(counts),
  gene_name = rownames(counts),
  feature_type = "Gene Expression",
  stringsAsFactors = FALSE
)
write.table(feat_df, file.path(out_dir, "features.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
cat("  → features.tsv\n")

# ── 4. Cell metadata ──
meta <- obj@meta.data
write.csv(meta, file.path(out_dir, "metadata.csv"), row.names = TRUE)
cat(sprintf("  → metadata.csv (%d cols × %d rows)\n", ncol(meta), nrow(meta)))

cat("\nDone! All files in:", out_dir, "\n")
cat("Ready for Python rebuilding.\n")
