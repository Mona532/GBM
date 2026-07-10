library(SeuratDisk)
library(Seurat)

ref <- readRDS("E:/GBM/results/GBM_annotation/annotated_seurat.rds")

# Merge fine cell types into broad categories (following paper's 12-class approach)
merge_map <- c(
  "B cell"                                      = "B cell",
  "T cell"                                      = "T cell",
  "naive T cell"                                = "T cell",
  "CD8-positive, alpha-beta cytotoxic T cell"   = "T cell",
  "natural killer cell"                         = "NK",
  "microglial cell"                             = "Microglia",
  "activated microglia"                         = "Microglia",
  "disease-associated microglia"                = "Microglia",
  "tumor-associated macrophage"                 = "Myeloid",
  "monocyte"                                    = "Myeloid",
  "classical monocyte"                          = "Myeloid",
  "neutrophil"                                  = "Neutrophil",
  "glioblastoma stem cell"                      = "Tumor",
  "oligodendrocyte"                             = "Oligodendrocyte",
  "astrocyte"                                   = "Astrocyte",
  "pericyte"                                    = "Pericyte",
  "innate lymphoid cell"                        = "ILC",
  "proliferating cell"                          = "Proliferating"
)

ref$cell_type <- merge_map[as.character(ref$llm_annotation_res0.5)]
ref <- ref[, !is.na(ref$cell_type)]  # remove Unknown
cat("Merged cell types:\n")
print(table(ref$cell_type))

# Keep only raw counts for cell2location
DefaultAssay(ref) <- "RNA"

dir.create("E:/GBM/results/cell2loc", showWarnings=FALSE, recursive=TRUE)
SaveH5Seurat(ref, "E:/GBM/results/cell2loc/ref_merged.h5Seurat", overwrite=TRUE)
Convert("E:/GBM/results/cell2loc/ref_merged.h5Seurat", dest="h5ad", overwrite=TRUE)
cat("ref_merged.h5ad saved\n")
