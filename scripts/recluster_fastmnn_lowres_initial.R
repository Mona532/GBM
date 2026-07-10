options(stringsAsFactors = FALSE)

.libPaths(c("E:/GBM/R/R-4.3.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

obj_path <- "E:/GBM/GBM_DATA/rerun_qc_fixed_fastmnn/GBM_4samples_qc_doublet_filtered_fastmnn.rds"
out_dir <- "E:/GBM/GBM_DATA/rerun_qc_fixed_fastmnn_lowres"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"

obj <- FindNeighbors(obj, reduction = "mnn", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.3, verbose = FALSE)

saveRDS(obj, file.path(out_dir, "GBM_4samples_qc_doublet_filtered_fastmnn_lowres.rds"))
write.csv(obj@meta.data, file.path(out_dir, "integrated_metadata_lowres.csv"), row.names = TRUE)

markers <- FindAllMarkers(
  obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)
write.csv(markers, file.path(out_dir, "cluster_markers_findall.csv"), row.names = FALSE)

top5 <- do.call(
  rbind,
  lapply(split(markers, markers$cluster), function(df) {
    df <- df[order(-df$avg_log2FC, -df$pct.1), ]
    head(df, 5)
  })
)
write.csv(top5, file.path(out_dir, "cluster_markers_top5.csv"), row.names = FALSE)

p_sample <- DimPlot(obj, reduction = "umap", group.by = "sample")
ggsave(file.path(out_dir, "umap_by_sample.png"), p_sample, width = 10, height = 7, dpi = 300)

p_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(file.path(out_dir, "umap_by_cluster.png"), p_cluster, width = 10, height = 7, dpi = 300)

cat("Low-resolution reclustering complete\n")
cat("Clusters:", length(unique(obj$seurat_clusters)), "\n")
