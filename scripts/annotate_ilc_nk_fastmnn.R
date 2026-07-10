options(stringsAsFactors = FALSE)

.libPaths(c("E:/GBM/R/R-4.3.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

obj_path <- "E:/GBM/GBM_DATA/integrated_fastmnn/GBM_4samples_integrated_fastmnn.rds"
out_dir <- "E:/GBM/GBM_DATA/integrated_fastmnn"

obj <- readRDS(obj_path)
clusters <- as.character(obj@meta.data[["seurat_clusters"]])

ilc_nk_map <- c(
  "5" = "ILC_NKT_like",
  "9" = "NK_T_like_cytotoxic",
  "13" = "NK_proliferating_like",
  "16" = "NK_cytotoxic",
  "17" = "ILC_NKT_like",
  "18" = "NK_cytotoxic",
  "21" = "ILC_like_KLRB1_IL7R",
  "22" = "ILC_NKT_like"
)

obj$ilc_nk_annotation <- "Other"
hit <- clusters %in% names(ilc_nk_map)
obj$ilc_nk_annotation[hit] <- unname(ilc_nk_map[clusters[hit]])
obj$ilc_nk_lineage <- ifelse(obj$ilc_nk_annotation == "Other", "Non_ILC_NK", "ILC_NK")

saveRDS(obj, file.path(out_dir, "GBM_4samples_integrated_fastmnn_annotated.rds"))

annot_counts <- as.data.frame(table(
  cluster = obj@meta.data[["seurat_clusters"]],
  ilc_nk_annotation = obj$ilc_nk_annotation
))
annot_counts <- annot_counts[annot_counts$Freq > 0, ]
write.csv(annot_counts, file.path(out_dir, "ilc_nk_annotation_counts.csv"), row.names = FALSE)

genes <- c("NKG7", "GNLY", "KLRD1", "KLRB1", "TRAC", "CD3D", "IL7R", "XCL1", "XCL2", "GZMB", "PRF1", "GATA3", "TRDC")

p1 <- DimPlot(obj, reduction = "umap", group.by = "ilc_nk_annotation", label = TRUE, repel = TRUE)
ggsave(
  file.path(out_dir, "umap_ilc_nk_annotation.png"),
  p1,
  width = 11,
  height = 8,
  dpi = 300
)

p2 <- DimPlot(obj, reduction = "umap", cells.highlight = WhichCells(obj, expression = ilc_nk_lineage == "ILC_NK"))
ggsave(
  file.path(out_dir, "umap_ilc_nk_highlight.png"),
  p2,
  width = 11,
  height = 8,
  dpi = 300
)

p3 <- DotPlot(obj, features = genes, group.by = "ilc_nk_annotation") + RotatedAxis()
ggsave(
  file.path(out_dir, "dotplot_ilc_nk_markers.png"),
  p3,
  width = 13,
  height = 7,
  dpi = 300
)

ilc_nk_obj <- subset(obj, subset = ilc_nk_lineage == "ILC_NK")
saveRDS(ilc_nk_obj, file.path(out_dir, "GBM_4samples_ILC_NK_subset.rds"))

cat("ILC/NK annotation complete\n")
cat("Annotated cells:", ncol(ilc_nk_obj), "\n")
