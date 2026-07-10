args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript extract_rds_to_10x.R <input.rds> <output_dir>")
}

input_rds <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input is not a Seurat object: ", input_rds)
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

barcodes <- colnames(counts)
if (is.null(barcodes)) {
  barcodes <- sprintf("cell_%d", seq_len(ncol(counts)))
  colnames(counts) <- barcodes
}

gene_names <- rownames(counts)
if (is.null(gene_names)) {
  gene_names <- sprintf("gene_%d", seq_len(nrow(counts)))
  rownames(counts) <- gene_names
}

meta <- obj@meta.data
meta <- meta[barcodes, , drop = FALSE]
meta$barcode <- rownames(meta)
meta <- meta[, c("barcode", setdiff(colnames(meta), "barcode")), drop = FALSE]

features <- data.frame(
  gene_id = gene_names,
  gene_name = gene_names,
  feature_type = rep("Gene Expression", length(gene_names)),
  stringsAsFactors = FALSE
)

matrix_path <- file.path(output_dir, "matrix.mtx")
features_path <- file.path(output_dir, "features.tsv")
barcodes_path <- file.path(output_dir, "barcodes.tsv")
metadata_path <- file.path(output_dir, "metadata.csv")

writeMM(counts, matrix_path)
write.table(
  features,
  file = features_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  data.frame(barcode = barcodes, stringsAsFactors = FALSE),
  file = barcodes_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.csv(meta, file = metadata_path, row.names = FALSE)

cat("Wrote:", output_dir, "\n")
cat("Cells:", ncol(counts), "Genes:", nrow(counts), "\n")
