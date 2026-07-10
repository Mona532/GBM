options(stringsAsFactors = FALSE)

.libPaths(c("E:/GBM/R/R-4.3.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

obj_path <- "E:/GBM/GBM_DATA/rerun_qc_fixed_fastmnn/GBM_4samples_qc_doublet_filtered_fastmnn.rds"
out_dir <- "E:/GBM/GBM_DATA/rerun_qc_fixed_fastmnn"

obj <- readRDS(obj_path)
clusters <- as.character(obj@meta.data[["seurat_clusters"]])

broad_map <- c(
  "0" = "SPP1_Macrophage",
  "1" = "Myeloid",
  "2" = "Microglia_like",
  "3" = "Myeloid",
  "4" = "Neutrophil",
  "5" = "NK_ILC",
  "6" = "Tumor_like_Neural",
  "7" = "Myeloid",
  "8" = "Tumor_like_Neural",
  "9" = "Oligodendrocyte",
  "10" = "T_cell",
  "11" = "Neutrophil",
  "12" = "Myeloid",
  "13" = "NK_ILC",
  "14" = "Myeloid",
  "15" = "Myeloid",
  "16" = "NK_ILC",
  "17" = "T_cell",
  "18" = "NK_ILC",
  "19" = "NK_ILC",
  "20" = "NK_ILC",
  "21" = "Myeloid",
  "22" = "Myeloid",
  "23" = "B_cell",
  "24" = "Oligodendrocyte"
)

obj$broad_celltype <- unname(broad_map[clusters])
obj$broad_celltype[is.na(obj$broad_celltype)] <- "Unassigned"
obj$nk_ilc_flag <- ifelse(obj$broad_celltype == "NK_ILC", "NK_ILC", "Non_NK_ILC")

saveRDS(obj, file.path(out_dir, "GBM_4samples_qc_doublet_filtered_fastmnn_annotated.rds"))

annot_counts <- as.data.frame(table(
  cluster = obj@meta.data[["seurat_clusters"]],
  broad_celltype = obj$broad_celltype
))
annot_counts <- annot_counts[annot_counts$Freq > 0, ]
write.csv(annot_counts, file.path(out_dir, "broad_celltype_counts.csv"), row.names = FALSE)

p1 <- DimPlot(obj, reduction = "umap", group.by = "broad_celltype", label = TRUE, repel = TRUE)
ggsave(
  file.path(out_dir, "umap_broad_celltype.png"),
  p1,
  width = 12,
  height = 8,
  dpi = 300
)

p2 <- DimPlot(obj, reduction = "umap", cells.highlight = WhichCells(obj, expression = nk_ilc_flag == "NK_ILC"))
ggsave(
  file.path(out_dir, "umap_nk_ilc_highlight.png"),
  p2,
  width = 12,
  height = 8,
  dpi = 300
)

genes <- c(
  "NKG7", "GNLY", "KLRD1", "KLRB1", "TRAC", "CD3D", "IL7R",
  "LYZ", "FCER1G", "TYROBP", "SPP1", "HLA-DRA", "FCN1",
  "MS4A1", "CD79A", "MRC1", "CD163", "PLP1", "MBP", "PTPRZ1"
)
p3 <- DotPlot(obj, features = genes, group.by = "broad_celltype") + RotatedAxis()
ggsave(
  file.path(out_dir, "dotplot_broad_celltype_markers.png"),
  p3,
  width = 13,
  height = 7,
  dpi = 300
)

cat("Broad annotation complete\n")
cat("NK/ILC cells:", sum(obj$broad_celltype == "NK_ILC"), "\n")
