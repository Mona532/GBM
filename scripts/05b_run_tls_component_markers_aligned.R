library(edgeR)

root <- "E:/GBM/results"
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
comp <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
weights <- read.csv(file.path(root, "tls_compnmf_rank4_unit_weights.csv"), check.names = FALSE)
weights$dominant_ecotype <- factor(weights$dominant_ecotype, levels = c("E1", "E2", "E3", "E4"))

if (is.null(rownames(comp))) {
  stop("Component matrix is missing rownames; expected unit_id rownames from pseudobulk build step.")
}

weights <- weights[match(rownames(comp), weights$unit_id), , drop = FALSE]
if (any(is.na(weights$unit_id))) {
  stop("Failed to align component weights to component matrix rownames.")
}

comp_df <- as.data.frame(comp)
comp_df$unit_id <- rownames(comp)

ann <- merge(
  weights[, c("unit_id", "sample", "n_spots", "tls_score_mean", "dominant_ecotype", "dominant_weight")],
  comp_df,
  by = "unit_id"
)

ann$maturity_score <- ann$B + ann$Plasma + ann$CD4_T + ann$CD8_T + ann$Dendritic + ann$NK
ann$ILC_total <- ann$ILC1 + ann$ILC2 + ann$ILC3
ann$ILC1_frac <- ann$ILC1 / pmax(ann$ILC_total, 1e-8)
ann$ILC2_frac <- ann$ILC2 / pmax(ann$ILC_total, 1e-8)
ann$ILC3_frac <- ann$ILC3 / pmax(ann$ILC_total, 1e-8)

celltypes <- c(
  "B", "CD4_T", "CD8_T", "Dendritic", "Glial", "Glioma",
  "ILC1", "ILC2", "ILC3", "Macrophage", "NK", "Plasma", "Vascular"
)

ec_sum <- do.call(rbind, lapply(levels(ann$dominant_ecotype), function(e) {
  sub <- ann[ann$dominant_ecotype == e, , drop = FALSE]
  means <- colMeans(sub[, celltypes, drop = FALSE])
  ord <- order(means, decreasing = TRUE)
  data.frame(
    ecotype = e,
    n_units = nrow(sub),
    n_samples = length(unique(sub$sample)),
    median_spots = median(sub$n_spots),
    mean_tls_score = mean(sub$tls_score_mean),
    mean_dominant_weight = mean(sub$dominant_weight),
    maturity_score = mean(sub$maturity_score),
    ILC_total = mean(sub$ILC_total),
    ILC1_frac = median(sub$ILC1_frac),
    ILC2_frac = median(sub$ILC2_frac),
    ILC3_frac = median(sub$ILC3_frac),
    top1 = names(means)[ord[1]],
    top2 = names(means)[ord[2]],
    top3 = names(means)[ord[3]],
    B = means["B"],
    CD4_T = means["CD4_T"],
    CD8_T = means["CD8_T"],
    Dendritic = means["Dendritic"],
    Glial = means["Glial"],
    Glioma = means["Glioma"],
    ILC1 = means["ILC1"],
    ILC2 = means["ILC2"],
    ILC3 = means["ILC3"],
    Macrophage = means["Macrophage"],
    NK = means["NK"],
    Plasma = means["Plasma"],
    Vascular = means["Vascular"],
    stringsAsFactors = FALSE
  )
}))
write.csv(ec_sum, file.path(root, "tls_compnmf_rank4_ecotype_annotated_summary.csv"), row.names = FALSE)

count_mat <- counts[, weights$unit_id, drop = FALSE]
group <- weights$dominant_ecotype
keep <- filterByExpr(count_mat, group = group)
y <- DGEList(counts = count_mat[keep, , drop = FALSE], group = group)
y <- calcNormFactors(y)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

for (e in levels(group)) {
  others <- setdiff(levels(group), e)
  coef_vec <- rep(-1 / length(others), length(levels(group)))
  names(coef_vec) <- levels(group)
  coef_vec[e] <- 1
  qlf <- glmQLFTest(fit, contrast = coef_vec)
  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$gene <- rownames(tab)
  write.csv(tab, file.path(root, paste0("tls_compnmf_rank4_", e, "_markers_vs_rest.csv")), row.names = FALSE)

  top_up <- tab[tab$logFC > 0, c("gene", "logFC", "FDR")]
  write.csv(head(top_up, 30), file.path(root, paste0("tls_compnmf_rank4_", e, "_top30_up.csv")), row.names = FALSE)
}

write.csv(
  ann[, c("unit_id", "sample", "dominant_ecotype", "maturity_score", "ILC_total", "ILC1", "ILC2", "ILC3")],
  file.path(root, "tls_component_features_per_component.csv"),
  row.names = FALSE
)

cat("Saved rank-4 TLS ecotype summaries and marker tables\n")
print(ec_sum)
