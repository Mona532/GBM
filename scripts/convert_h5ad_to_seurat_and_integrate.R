args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript convert_h5ad_to_seurat_and_integrate.R <input_h5ad> <output_rds> <output_integrated_rds>")
}

input_h5ad <- normalizePath(args[[1]], mustWork = TRUE)
output_rds <- args[[2]]
output_integrated_rds <- args[[3]]

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
  library(Seurat)
})

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_integrated_rds), recursive = TRUE, showWarnings = FALSE)

sce <- readH5AD(input_h5ad, use_hdf5 = FALSE)

counts <- assay(sce, "X")
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

meta <- as.data.frame(colData(sce))
meta$cell_id <- colnames(sce)
rownames(meta) <- meta$cell_id

if (!"sample" %in% colnames(meta)) {
  stop("Missing required obs column: sample")
}

if (!"Source" %in% colnames(meta)) {
  meta$Source <- ""
}

features <- data.frame(
  gene_symbol = rownames(sce),
  stringsAsFactors = FALSE,
  row.names = rownames(sce)
)

obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  assay = "RNA",
  project = "GBM"
)
obj[["RNA"]]@meta.features <- features

saveRDS(obj, file = output_rds)

objs <- SplitObject(obj, split.by = "sample")
objs <- lapply(
  objs,
  function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    x <- ScaleData(x, verbose = FALSE)
    x <- RunPCA(x, npcs = 30, verbose = FALSE)
    x
  }
)

integration_features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
anchors <- FindIntegrationAnchors(
  object.list = objs,
  anchor.features = integration_features,
  reduction = "rpca"
)
integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 1.0, verbose = FALSE)

DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)

saveRDS(integrated, file = output_integrated_rds)

cat("Wrote raw Seurat RDS:", output_rds, "\n")
cat("Wrote integrated Seurat RDS:", output_integrated_rds, "\n")
cat("Cells:", ncol(integrated), "Genes:", nrow(integrated), "\n")
