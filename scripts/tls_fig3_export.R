# Export per-ecotype logCPM gene expression for Fig3 heatmap
root <- "E:/GBM/results"
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
meta <- read.csv(file.path(root, "tls_pseudobulk_component_metadata.csv"))
nmf <- read.csv(file.path(root, "tls_compnmf_rank5_unit_weights.csv"), check.names=FALSE)
markers <- read.csv("E:/GBM/Pan-Cancer_Spatial_Atlas_TLS-main/data/fig3_markers.csv")

eco_names <- c("E1"="Glial-CD4","E2"="TLS-structural","E3"="Vascular","E4"="Lymphocyte","E5"="Myeloid")
eco_order <- c("Lymphocyte","TLS-structural","Glial-CD4","Vascular","Myeloid")

common <- intersect(colnames(counts), nmf$unit_id)
counts <- counts[, common, drop=FALSE]
comp_eco <- nmf$dominant_ecotype[match(common, nmf$unit_id)]
comp_eco_name <- eco_names[comp_eco]

# CPM + log2
X <- as.matrix(counts)
X_cpm <- sweep(X, 2, colSums(X), "/") * 1e6
X_logcpm <- log2(X_cpm + 1)

avail <- intersect(markers$gene, rownames(counts))
cat(sprintf("Genes: %d/%d\n", length(avail), nrow(markers)))

eco_means <- matrix(NA, length(eco_order), length(avail))
rownames(eco_means) <- eco_order
colnames(eco_means) <- avail
for (i in seq_along(eco_order)) {
  mask <- comp_eco_name == eco_order[i]
  eco_means[i, ] <- rowMeans(X_logcpm[avail, mask, drop=FALSE])
}

eco_scaled <- scale(eco_means)
eco_scaled <- t(eco_scaled)

saveRDS(eco_scaled, file.path(root, "tls_fig3_expr_scaled.rds"))
write.csv(eco_scaled, file.path(root, "tls_fig3_expr_scaled.csv"))
cat(sprintf("Saved: %d genes x %d ecotypes\n", nrow(eco_scaled), ncol(eco_scaled)))
