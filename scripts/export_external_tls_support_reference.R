suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

out_dir <- 'E:/GBM/results/reference_rebuild/shareable_unified_reference/external_tls_support_export'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

get_counts <- function(obj) {
  if ('RNA' %in% names(obj@assays)) {
    tryCatch(GetAssayData(obj, assay = 'RNA', slot = 'counts'), error = function(e) GetAssayData(obj, assay = 'RNA', layer = 'counts'))
  } else {
    stop('RNA assay not found')
  }
}

sample_cells <- function(cells, cap) {
  if (length(cells) <= cap) return(cells)
  sample(cells, cap)
}

parts <- list()
meta_parts <- list()

obj_t <- readRDS('E:/GBM/results/external_sc_integration/tcell_recluster/tcell_lineage_harmony_clustered.rds')
md_t <- obj_t@meta.data
cells_treg <- rownames(md_t)[as.character(md_t$seurat_clusters) == '4']
cells_treg <- sample_cells(cells_treg, 250)
counts_t <- get_counts(obj_t)[, cells_treg, drop = FALSE]
meta_t <- md_t[cells_treg, c('sample_id','dataset_id','Cancer.type','Tissue','Celltype'), drop = FALSE]
meta_t$ref_label <- 'Treg'
meta_t$ref_source <- 'external_tcell_recluster'
parts[['Treg']] <- counts_t
meta_parts[['Treg']] <- meta_t

obj_13 <- readRDS('E:/GBM/results/external_sc_integration/targeted_recluster/cluster_13/cluster_13_reclustered.rds')
ann_13 <- fread('E:/GBM/results/external_sc_integration/targeted_recluster/cluster_13/cluster_13_provisional_annotation.csv')
map_13 <- setNames(ann_13$label, as.character(ann_13$subcluster))
md_13 <- obj_13@meta.data
lab_13 <- map_13[as.character(md_13$seurat_clusters)]
cells_tfh <- rownames(md_13)[lab_13 == 'activated_Tfh_candidate']
counts_13 <- get_counts(obj_13)[, cells_tfh, drop = FALSE]
meta_13 <- md_13[cells_tfh, c('sample_id','dataset_id','Cancer.type','Tissue','Celltype'), drop = FALSE]
meta_13$ref_label <- 'Tfh-like_CD4'
meta_13$ref_source <- 'external_cluster13_recluster'
parts[['Tfh-like_CD4']] <- counts_13
meta_parts[['Tfh-like_CD4']] <- meta_13

obj_10 <- readRDS('E:/GBM/results/external_sc_integration/targeted_recluster/cluster_10/cluster_10_reclustered.rds')
ann_10 <- fread('E:/GBM/results/external_sc_integration/targeted_recluster/cluster_10/cluster_10_provisional_annotation.csv')
map_10 <- setNames(ann_10$label, as.character(ann_10$subcluster))
md_10 <- obj_10@meta.data
lab_10 <- map_10[as.character(md_10$seurat_clusters)]
keep_dc <- lab_10 %in% c('cDC1_like', 'cDC2_migratory_like', 'DC_like_minor')
cells_dc <- rownames(md_10)[keep_dc]
cells_dc <- sample_cells(cells_dc, 300)
counts_10 <- get_counts(obj_10)[, cells_dc, drop = FALSE]
meta_10 <- md_10[cells_dc, c('sample_id','dataset_id','Cancer.type','Tissue','Celltype'), drop = FALSE]
meta_10$ref_label <- 'cDC_or_mature_DC'
meta_10$ref_source <- 'external_cluster10_recluster'
parts[['cDC_or_mature_DC']] <- counts_10
meta_parts[['cDC_or_mature_DC']] <- meta_10

obj_s <- readRDS('E:/GBM/results/external_sc_integration/stromal_endothelial_recluster/stromal_endothelial_harmony_clustered.rds')
ann_s <- fread('E:/GBM/results/external_sc_integration/stromal_endothelial_recluster/provisional_stromal_endothelial_annotation.csv')
hev_clusters <- ann_s$cluster[ann_s$label == 'HEV_like_endothelial_candidate']
md_s <- obj_s@meta.data
cells_hev <- rownames(md_s)[as.character(md_s$seurat_clusters) %in% as.character(hev_clusters)]
cells_hev <- sample_cells(cells_hev, 150)
counts_s <- get_counts(obj_s)[, cells_hev, drop = FALSE]
meta_s <- md_s[cells_hev, c('sample_id','dataset_id','Cancer.type','Tissue','Celltype'), drop = FALSE]
meta_s$ref_label <- 'HEV-like_endothelial'
meta_s$ref_source <- 'external_stromal_endothelial_recluster'
parts[['HEV-like_endothelial']] <- counts_s
meta_parts[['HEV-like_endothelial']] <- meta_s

common_genes <- Reduce(intersect, lapply(parts, rownames))
parts <- lapply(parts, function(x) x[common_genes, , drop = FALSE])
combined <- do.call(cbind, parts)
meta <- rbindlist(meta_parts, use.names = TRUE, fill = TRUE, idcol = 'part_name')
meta$cell_barcode <- colnames(combined)
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_barcode
meta <- meta[colnames(combined), , drop = FALSE]

writeMM(combined, file.path(out_dir, 'matrix.mtx'))
write.table(data.frame(gene = common_genes), file.path(out_dir, 'features.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(cell = colnames(combined)), file.path(out_dir, 'barcodes.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(meta, file.path(out_dir, 'metadata.csv'), row.names = FALSE)
write.csv(as.data.frame(table(meta$ref_label)), file.path(out_dir, 'label_counts.csv'), row.names = FALSE)
cat('saved', out_dir, '\n')
cat('dims', nrow(combined), ncol(combined), '\n')
print(table(meta$ref_label))
