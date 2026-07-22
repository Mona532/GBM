options(stringsAsFactors = FALSE)

.libPaths(c("E:/GBM/R/R-4.3.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(batchelor)
  library(ggplot2)
})

data_dir <- "E:/GBM/GBM_DATA"
out_dir <- file.path(data_dir, "rerun_qc_fixed_fastmnn")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

apply_fixed_qc_and_doublet <- function(obj, sample_id) {
  DefaultAssay(obj) <- "RNA"
  obj$sample <- sample_id
  obj$dataset <- sample_id
  obj$batch <- sample_id
  obj$source_dataset <- sample_id

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]")

  md <- obj@meta.data
  keep <- md[["nFeature_RNA"]] >= 200 &
    md[["nFeature_RNA"]] <= 7500 &
    md[["nCount_RNA"]] >= 500 &
    md[["nCount_RNA"]] <= 50000 &
    md[["percent.mt"]] < 20

  before_n <- ncol(obj)
  obj <- subset(obj, cells = colnames(obj)[keep])
  after_qc_n <- ncol(obj)

  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  set.seed(1234)
  sce <- scDblFinder(sce, samples = rep(sample_id, ncol(sce)))
  dbl_df <- as.data.frame(colData(sce))

  obj[["scDblFinder.sample"]] <- dbl_df[["scDblFinder.sample"]]
  obj[["scDblFinder.class"]] <- dbl_df[["scDblFinder.class"]]
  obj[["scDblFinder.score"]] <- dbl_df[["scDblFinder.score"]]
  obj[["scDblFinder.weighted"]] <- dbl_df[["scDblFinder.weighted"]]

  singlet_cells <- rownames(obj@meta.data)[obj@meta.data[["scDblFinder.class"]] == "singlet"]
  obj <- subset(obj, cells = singlet_cells)
  obj <- RenameCells(obj, add.cell.id = sample_id)

  summary_df <- data.frame(
    sample = sample_id,
    before_cells = before_n,
    after_qc_cells = after_qc_n,
    singlet_cells = ncol(obj),
    removed_qc = before_n - after_qc_n,
    removed_doublets = after_qc_n - ncol(obj),
    stringsAsFactors = FALSE
  )

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

  list(object = obj, summary = summary_df)
}

gbm01 <- readRDS(file.path(data_dir, "GBM01.rds"))
gbm02 <- readRDS(file.path(data_dir, "GBM02.rds"))

ilc1_10x <- Read10X(file.path(data_dir, "GBM-ILC1"))
ilc2_10x <- Read10X(file.path(data_dir, "GBM-ILC2"))

ilc1 <- CreateSeuratObject(ilc1_10x[["Gene Expression"]], assay = "RNA", project = "GBM-ILC1")
ilc2 <- CreateSeuratObject(ilc2_10x[["Gene Expression"]], assay = "RNA", project = "GBM-ILC2")

res01 <- apply_fixed_qc_and_doublet(gbm01, "GBM01")
res02 <- apply_fixed_qc_and_doublet(gbm02, "GBM02")
res11 <- apply_fixed_qc_and_doublet(ilc1, "GBM-ILC1")
res12 <- apply_fixed_qc_and_doublet(ilc2, "GBM-ILC2")

objs <- list(res01$object, res02$object, res11$object, res12$object)
qc_summary <- rbind(res01$summary, res02$summary, res11$summary, res12$summary)
write.csv(qc_summary, file.path(out_dir, "qc_doublet_summary.csv"), row.names = FALSE)

saveRDS(res01$object, file.path(out_dir, "GBM01_qc_singlets.rds"))
saveRDS(res02$object, file.path(out_dir, "GBM02_qc_singlets.rds"))
saveRDS(res11$object, file.path(out_dir, "GBM-ILC1_qc_singlets.rds"))
saveRDS(res12$object, file.path(out_dir, "GBM-ILC2_qc_singlets.rds"))

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
sce_list <- lapply(objs, function(obj) as.SingleCellExperiment(obj, assay = "RNA"))
mnn_out <- do.call(
  batchelor::fastMNN,
  c(sce_list, list(subset.row = features, d = 50, auto.merge = TRUE))
)

merged <- Reduce(function(x, y) merge(x, y), objs)
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
integrated <- FindClusters(integrated, resolution = 0.3, verbose = FALSE)

saveRDS(integrated, file.path(out_dir, "GBM_4samples_qc_doublet_filtered_fastmnn.rds"))
write.csv(integrated@meta.data, file.path(out_dir, "integrated_metadata.csv"), row.names = TRUE)

p_sample <- DimPlot(integrated, reduction = "umap", group.by = "sample")
ggsave(file.path(out_dir, "umap_by_sample.png"), p_sample, width = 10, height = 7, dpi = 300)

p_cluster <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(file.path(out_dir, "umap_by_cluster.png"), p_cluster, width = 10, height = 7, dpi = 300)

cat("QC + doublet removal + FastMNN complete\n")
cat("Output directory:", out_dir, "\n")
cat("Integrated cells:", ncol(integrated), "\n")
cat("Integrated genes:", nrow(integrated), "\n")
