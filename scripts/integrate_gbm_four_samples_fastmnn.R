options(stringsAsFactors = FALSE)

.libPaths(c("E:/GBM/R/R-4.3.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(batchelor)
  library(SingleCellExperiment)
  library(ggplot2)
})

data_dir <- "E:/GBM/GBM_DATA"
out_dir <- file.path(data_dir, "integrated_fastmnn")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

add_missing_meta_columns <- function(obj, required_cols) {
  for (col in required_cols) {
    if (!col %in% colnames(obj@meta.data)) {
      obj[[col]] <- NA
    }
  }
  obj
}

standardize_object <- function(obj, sample_id, dataset_label) {
  DefaultAssay(obj) <- "RNA"
  obj$sample <- sample_id
  obj$dataset <- dataset_label
  obj$batch <- sample_id
  obj$source_dataset <- dataset_label
  obj <- RenameCells(obj, add.cell.id = sample_id)
  obj
}

load_rds_object <- function(path, sample_id) {
  obj <- readRDS(path)
  standardize_object(obj, sample_id = sample_id, dataset_label = sample_id)
}

load_10x_object <- function(dir_path, sample_id) {
  x <- Read10X(dir_path)
  obj <- CreateSeuratObject(
    counts = x[["Gene Expression"]],
    assay = "RNA",
    project = sample_id
  )

  if ("Antibody Capture" %in% names(x)) {
    obj[["HTO"]] <- CreateAssayObject(counts = x[["Antibody Capture"]])
  }

  meta_path <- file.path(dir_path, "metadata.csv")
  if (file.exists(meta_path)) {
    meta <- read.csv(meta_path, check.names = FALSE)
    key_col <- if ("barcode" %in% colnames(meta)) "barcode" else colnames(meta)[1]
    rownames(meta) <- meta[[key_col]]
    meta <- meta[colnames(obj), setdiff(colnames(meta), key_col), drop = FALSE]
    obj <- AddMetaData(obj, metadata = meta)
  }

  standardize_object(obj, sample_id = sample_id, dataset_label = sample_id)
}

required_meta <- c(
  "orig.ident",
  "nCount_RNA",
  "nFeature_RNA",
  "nCount_HTO",
  "nFeature_HTO",
  "HTO_maxID",
  "HTO_secondID",
  "HTO_margin",
  "HTO_classification",
  "HTO_classification.global",
  "hash.ID",
  "RNA_snn_res.0.6",
  "seurat_clusters",
  "Source",
  "barcode_raw",
  "barcode_merged",
  "dataset",
  "sample",
  "batch",
  "source_dataset"
)

objs <- list(
  GBM01 = load_rds_object(file.path(data_dir, "GBM01.rds"), "GBM01"),
  GBM02 = load_rds_object(file.path(data_dir, "GBM02.rds"), "GBM02"),
  `GBM-ILC1` = load_10x_object(file.path(data_dir, "GBM-ILC1"), "GBM-ILC1"),
  `GBM-ILC2` = load_10x_object(file.path(data_dir, "GBM-ILC2"), "GBM-ILC2")
)

objs <- lapply(objs, add_missing_meta_columns, required_cols = required_meta)

sample_summary <- do.call(
  rbind,
  lapply(objs, function(obj) {
    data.frame(
      sample = unique(obj$sample)[1],
      cells = ncol(obj),
      genes = nrow(obj),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(sample_summary, file.path(out_dir, "sample_summary.csv"), row.names = FALSE)

objs <- lapply(
  objs,
  function(obj) {
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    obj
  }
)

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)

sce_list <- lapply(
  objs,
  function(obj) as.SingleCellExperiment(obj, assay = "RNA")
)

mnn_args <- c(
  sce_list,
  list(subset.row = features, d = 50, auto.merge = TRUE)
)
mnn_out <- do.call(batchelor::fastMNN, mnn_args)

merged <- Reduce(
  function(x, y) merge(x, y),
  objs
)
DefaultAssay(merged) <- "RNA"

corrected <- reducedDim(mnn_out, "corrected")
corrected <- corrected[colnames(merged), 1:min(50, ncol(corrected)), drop = FALSE]
merged[["mnn"]] <- CreateDimReducObject(
  embeddings = corrected,
  key = "MNN_",
  assay = "RNA"
)

integrated <- merged
integrated <- RunUMAP(integrated, reduction = "mnn", dims = 1:30, reduction.name = "umap", verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "mnn", dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.6, verbose = FALSE)

saveRDS(integrated, file.path(out_dir, "GBM_4samples_integrated_fastmnn.rds"))

umap_by_sample <- DimPlot(integrated, reduction = "umap", group.by = "sample")
ggsave(
  filename = file.path(out_dir, "GBM_4samples_fastmnn_umap_by_sample.png"),
  plot = umap_by_sample,
  width = 10,
  height = 7,
  dpi = 300
)

umap_by_cluster <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(
  filename = file.path(out_dir, "GBM_4samples_fastmnn_umap_by_cluster.png"),
  plot = umap_by_cluster,
  width = 10,
  height = 7,
  dpi = 300
)

write.csv(
  integrated@meta.data,
  file.path(out_dir, "GBM_4samples_fastmnn_metadata.csv"),
  row.names = TRUE
)

cat("FastMNN integration complete\n")
cat("Output directory:", out_dir, "\n")
cat("Integrated object cells:", ncol(integrated), "\n")
cat("Integrated object genes:", nrow(integrated), "\n")
