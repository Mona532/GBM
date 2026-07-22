suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(readxl)
})

root <- "E:/GBM/results"
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
weights <- read.csv(file.path(root, "tls_compnmf_rank5_unit_weights.csv"), check.names = FALSE, stringsAsFactors = FALSE)
comp <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
maturity_genes <- read.csv(file.path(root, "maturity_genes_used.csv"), check.names = FALSE, stringsAsFactors = FALSE)
eco_summary <- read.csv(file.path(root, "tls_compnmf_rank5_ecotype_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)

weights <- weights[match(colnames(counts), weights$unit_id), , drop = FALSE]
if (any(is.na(weights$unit_id))) {
  stop("Failed to align rank5 unit weights to count matrix.")
}
if (!identical(weights$unit_id, colnames(counts))) {
  stop("Weight order mismatch after alignment.")
}

weights$sample <- as.character(weights$sample)
weights$ecotype <- factor(weights$dominant_ecotype, levels = paste0("E", 1:5))
comp <- comp[weights$unit_id, , drop = FALSE]

rx <- read_excel("E:/GBM/ti tianran2.xlsx")
cat_names <- c("Glutamate", "GABA/Gly", "Cholinergic", "DA/NE", "Serotonin")
receptor_tbl <- do.call(rbind, lapply(seq_along(rx), function(i) {
  genes <- unique(na.omit(as.character(rx[[i]])))
  data.frame(gene = genes, category = cat_names[i], stringsAsFactors = FALSE)
}))
receptor_tbl <- receptor_tbl[!duplicated(receptor_tbl$gene), , drop = FALSE]
bad_marker_pattern <- "^(RPS|RPL|MT-|MTRNR|LINC|MALAT1$|NEAT1$)"

make_contrast <- function(target, levels_vec, coef_names, prefix) {
  v <- rep(0, length(coef_names))
  names(v) <- coef_names
  others <- setdiff(levels_vec, target)
  if (target == levels_vec[1]) {
    for (other in others) v[paste0(prefix, other)] <- -1 / length(others)
  } else {
    v[paste0(prefix, target)] <- 1
    nonbaseline_others <- setdiff(others, levels_vec[1])
    for (other in nonbaseline_others) v[paste0(prefix, other)] <- -1 / length(others)
  }
  v
}

zmat <- function(x) {
  x <- as.matrix(x)
  out <- t(scale(t(x)))
  out[!is.finite(out)] <- 0
  out
}

plot_heatmap <- function(mat, title, fill_name, out_jpg, out_pdf, low = "#2166ac", mid = "#f7f7f7", high = "#b2182b") {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("feature", "group", "value")
  df$feature <- factor(df$feature, levels = rev(rownames(mat)))
  df$group <- factor(df$group, levels = colnames(mat))
  p <- ggplot(df, aes(x = group, y = feature, fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0, name = fill_name) +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
      legend.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black")
    )
  ggsave(out_jpg, p, width = 8.5, height = max(4.5, nrow(mat) * 0.28), dpi = 300)
  ggsave(out_pdf, p, width = 8.5, height = max(4.5, nrow(mat) * 0.28))
}

plot_heatmap_continuous <- function(mat, title, fill_name, out_jpg, out_pdf, low = "#F6E8E3", high = "red4") {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("feature", "group", "value")
  df$feature <- factor(df$feature, levels = rev(rownames(mat)))
  df$group <- factor(df$group, levels = colnames(mat))
  p <- ggplot(df, aes(x = group, y = feature, fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    scale_fill_gradient(low = low, high = high, name = fill_name) +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
      legend.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black")
    )
  ggsave(out_jpg, p, width = 8.5, height = max(4.5, nrow(mat) * 0.28), dpi = 300)
  ggsave(out_pdf, p, width = 8.5, height = max(4.5, nrow(mat) * 0.28))
}

## 1) Sample-aware DEG
sample_eco_key <- paste(weights$sample, weights$ecotype, sep = "__")
split_idx <- split(seq_len(ncol(counts)), sample_eco_key)
agg_counts <- do.call(cbind, lapply(split_idx, function(idx) rowSums(counts[, idx, drop = FALSE])))
agg_meta <- do.call(rbind, lapply(names(split_idx), function(key) {
  idx <- split_idx[[key]]
  data.frame(
    agg_id = key,
    sample = weights$sample[idx[1]],
    ecotype = as.character(weights$ecotype[idx[1]]),
    n_components = length(idx),
    mean_dominant_weight = mean(weights$dominant_weight[idx]),
    stringsAsFactors = FALSE
  )
}))
agg_meta$ecotype <- factor(agg_meta$ecotype, levels = levels(weights$ecotype))
colnames(agg_counts) <- agg_meta$agg_id

keep <- filterByExpr(agg_counts, group = agg_meta$ecotype)
y <- DGEList(counts = agg_counts[keep, , drop = FALSE], group = agg_meta$ecotype)
y <- calcNormFactors(y)
design <- model.matrix(~ sample + ecotype, data = agg_meta)
v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)

deg_top_genes <- list()
for (eco in levels(agg_meta$ecotype)) {
  contrast <- make_contrast(eco, levels(agg_meta$ecotype), colnames(design), "ecotype")
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  tab <- topTable(fit2, number = Inf, sort.by = "P")
  tab$gene <- rownames(tab)
  tab$ecotype <- eco
  tab <- tab[, c("ecotype", "gene", setdiff(names(tab), c("ecotype", "gene")))]
  write.csv(tab, file.path(root, paste0("tls_rank5_sample_aware_", eco, "_markers_vs_rest.csv")), row.names = FALSE)

  top_up <- tab[tab$logFC > 0, c("ecotype", "gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")]
  write.csv(head(top_up, 30), file.path(root, paste0("tls_rank5_sample_aware_", eco, "_top30_up.csv")), row.names = FALSE)

  plot_top <- top_up[!grepl(bad_marker_pattern, top_up$gene, ignore.case = TRUE), , drop = FALSE]
  sig_top <- plot_top[plot_top$adj.P.Val < 0.1, , drop = FALSE]
  if (nrow(sig_top) == 0) sig_top <- plot_top
  if (nrow(sig_top) == 0) sig_top <- top_up
  deg_top_genes[[eco]] <- head(sig_top$gene, 8)
}

logcpm_agg <- cpm(y, log = TRUE, prior.count = 2)
eco_means <- sapply(levels(agg_meta$ecotype), function(eco) {
  cols <- agg_meta$agg_id[agg_meta$ecotype == eco]
  rowMeans(logcpm_agg[, cols, drop = FALSE])
})
top_gene_vec <- unique(unlist(deg_top_genes, use.names = FALSE))
top_gene_vec <- top_gene_vec[!grepl(bad_marker_pattern, top_gene_vec, ignore.case = TRUE)]
top_gene_vec <- top_gene_vec[top_gene_vec %in% rownames(eco_means)]
deg_heat <- eco_means[top_gene_vec, , drop = FALSE]
deg_heat_z <- zmat(deg_heat)
write.csv(cbind(gene = rownames(deg_heat_z), deg_heat_z), file.path(root, "tls_rank5_deg_heatmap_matrix.csv"), row.names = FALSE)
plot_heatmap(
  deg_heat_z,
  "Rank5 ecotype DEG heatmap",
  "Row z-score",
  file.path(root, "fig_tls_rank5_deg_heatmap.jpg"),
  file.path(root, "fig_tls_rank5_deg_heatmap.pdf"),
  low = "#559B80",
  mid = "white",
  high = "red4"
)

## Common matrices
logcpm_comp <- cpm(DGEList(counts = counts), log = TRUE, prior.count = 2)
cpm_comp <- cpm(DGEList(counts = counts), log = FALSE)

## 2) E4 ILC1/2/3 vs receptor correlation heatmap
e4_units <- weights$unit_id[weights$ecotype == "E4"]
if (length(e4_units) < 6) stop("Too few E4 components for correlation analysis.")
ilc_cols <- c("ILC1", "ILC2", "ILC3")
ilc_raw <- as.matrix(comp[e4_units, ilc_cols, drop = FALSE])
cor_rows <- list()
for (gene in receptor_tbl$gene) {
  if (!gene %in% rownames(logcpm_comp)) next
  expr <- as.numeric(logcpm_comp[gene, e4_units])
  if (sd(expr) == 0) next
  for (ilc in ilc_cols) {
    vals <- as.numeric(ilc_raw[, ilc])
    if (sd(vals) == 0) next
    ct <- suppressWarnings(cor.test(expr, vals, method = "spearman", exact = FALSE))
    cor_rows[[length(cor_rows) + 1L]] <- data.frame(
      gene = gene,
      category = receptor_tbl$category[match(gene, receptor_tbl$gene)],
      ilc = ilc,
      rho = unname(ct$estimate),
      pvalue = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
}
cor_df <- do.call(rbind, cor_rows)
cor_df$fdr <- ave(cor_df$pvalue, cor_df$ilc, FUN = p.adjust, method = "BH")
write.csv(cor_df, file.path(root, "tls_rank5_E4_ilc_receptor_correlation.csv"), row.names = FALSE)

cor_df$category <- factor(cor_df$category, levels = cat_names)
for (cat in cat_names) {
  sub <- cor_df[cor_df$category == cat, , drop = FALSE]
  if (nrow(sub) == 0) next
  ord <- aggregate(abs(rho) ~ gene, data = sub, FUN = max)
  ord <- ord[order(ord$`abs(rho)`, decreasing = TRUE), "gene"]
  cor_mat <- xtabs(rho ~ gene + ilc, data = sub)
  cor_mat <- cor_mat[match(ord, rownames(cor_mat)), ilc_cols, drop = FALSE]
  out_stub <- gsub("[^A-Za-z0-9]+", "_", cat)
  write.csv(
    cbind(gene = rownames(cor_mat), cor_mat),
    file.path(root, paste0("tls_rank5_E4_ilc_receptor_heatmap_matrix_", out_stub, ".csv")),
    row.names = FALSE
  )
  plot_heatmap(
    cor_mat,
    paste0("E4 ", cat, " receptors vs ILC abundance"),
    "Spearman rho",
    file.path(root, paste0("fig_tls_rank5_E4_ilc_receptor_heatmap_", out_stub, ".jpg")),
    file.path(root, paste0("fig_tls_rank5_E4_ilc_receptor_heatmap_", out_stub, ".pdf")),
    low = "#559B80",
    mid = "white",
    high = "red4"
  )
}

## 3) Maturity heatmap
gene_sets <- split(maturity_genes$gene[maturity_genes$present == "yes"], maturity_genes$module[maturity_genes$present == "yes"])
module_mean <- function(genes) {
  genes <- intersect(genes, rownames(logcpm_comp))
  if (length(genes) == 0) return(rep(NA_real_, ncol(logcpm_comp)))
  colMeans(logcpm_comp[genes, , drop = FALSE])
}
score_df <- data.frame(
  unit_id = colnames(logcpm_comp),
  ecotype = weights$ecotype,
  B_diff_score = module_mean(gene_sets[["B diff (+)"]]) - module_mean(gene_sets[["B diff (-)"]]),
  Ig_production = module_mean(gene_sets[["Ig production"]]),
  GC_markers = module_mean(gene_sets[["GC markers"]]),
  Tfh_program = module_mean(gene_sets[["Tfh"]]),
  CXCL13_CCL19_CCL21 = module_mean(gene_sets[["CXCL chemokines"]]),
  B_abund = comp[colnames(logcpm_comp), "B"],
  Plasma_abund = comp[colnames(logcpm_comp), "Plasma"],
  Tfh_abund = comp[colnames(logcpm_comp), "Tfh-like_CD4"],
  FDC_abund = comp[colnames(logcpm_comp), "FDC"],
  HEV_abund = comp[colnames(logcpm_comp), "HEV-like_endothelial"],
  stringsAsFactors = FALSE
)

score_df$sample <- weights$sample[match(score_df$unit_id, weights$unit_id)]
score_df$sample_ecotype <- paste(score_df$sample, score_df$ecotype, sep = "__")
maturity_features <- setdiff(colnames(score_df), c("unit_id", "sample", "ecotype", "sample_ecotype"))

agg_score <- aggregate(score_df[, maturity_features, drop = FALSE], by = list(sample_ecotype = score_df$sample_ecotype), FUN = mean, na.rm = TRUE)
agg_score$sample <- sub("__.*$", "", agg_score$sample_ecotype)
agg_score$ecotype <- sub("^.*__", "", agg_score$sample_ecotype)
agg_score$ecotype <- factor(agg_score$ecotype, levels = levels(weights$ecotype))

maturity_mean <- sapply(levels(agg_score$ecotype), function(eco) {
  idx <- agg_score$ecotype == eco
  colMeans(agg_score[idx, maturity_features, drop = FALSE], na.rm = TRUE)
})

write.csv(cbind(feature = rownames(maturity_mean), maturity_mean), file.path(root, "tls_rank5_maturity_heatmap_raw_matrix.csv"), row.names = FALSE)
maturity_z <- zmat(maturity_mean)
write.csv(cbind(feature = rownames(maturity_z), maturity_z), file.path(root, "tls_rank5_maturity_heatmap_matrix.csv"), row.names = FALSE)
plot_heatmap(
  maturity_z,
  "Rank5 ecotype maturity score heatmap",
  "Row z-score",
  file.path(root, "fig_tls_rank5_maturity_heatmap.jpg"),
  file.path(root, "fig_tls_rank5_maturity_heatmap.pdf"),
  low = "#559B80",
  mid = "white",
  high = "red4"
)

## 4) Receptor dotplot across ecotypes
dot_rows <- list()
for (gene in receptor_tbl$gene) {
  if (!gene %in% rownames(logcpm_comp)) next
  for (eco in levels(weights$ecotype)) {
    units <- weights$unit_id[weights$ecotype == eco]
    if (length(units) == 0) next
    dot_rows[[length(dot_rows) + 1L]] <- data.frame(
      gene = gene,
      category = receptor_tbl$category[match(gene, receptor_tbl$gene)],
      ecotype = eco,
      detect_rate = mean(cpm_comp[gene, units] > 1),
      mean_expr = mean(logcpm_comp[gene, units]),
      stringsAsFactors = FALSE
    )
  }
}
dot_df <- do.call(rbind, dot_rows)
write.csv(dot_df, file.path(root, "tls_rank5_receptor_dotplot.csv"), row.names = FALSE)

gene_rank <- aggregate(mean_expr ~ gene + category, data = dot_df, FUN = max)
gene_rank <- gene_rank[order(gene_rank$category, gene_rank$mean_expr, decreasing = TRUE), ]
dot_df$ecotype <- factor(dot_df$ecotype, levels = levels(weights$ecotype))
for (cat in cat_names) {
  sub <- dot_df[dot_df$category == cat, , drop = FALSE]
  gene_order <- aggregate(mean_expr ~ gene, data = sub, FUN = max)
  gene_order <- gene_order[order(gene_order$mean_expr, decreasing = TRUE), "gene"]
  sub$gene <- factor(sub$gene, levels = rev(gene_order))
  out_stub <- gsub("[^A-Za-z0-9]+", "_", cat)
  p <- ggplot(sub, aes(x = ecotype, y = gene)) +
    geom_point(aes(size = detect_rate, color = mean_expr), alpha = 0.92) +
    scale_size(range = c(1.2, 7), name = "Detection") +
    scale_color_gradient(low = "#F6E8E3", high = "red4", name = "Mean logCPM") +
    labs(x = NULL, y = NULL, title = cat) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 7, colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
      legend.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black")
    )
  ggsave(
    file.path(root, paste0("fig_tls_rank5_receptor_dotplot_", out_stub, ".jpg")),
    p,
    width = 4.6,
    height = max(4.2, length(unique(sub$gene)) * 0.18),
    dpi = 300
  )
  ggsave(
    file.path(root, paste0("fig_tls_rank5_receptor_dotplot_", out_stub, ".pdf")),
    p,
    width = 4.6,
    height = max(4.2, length(unique(sub$gene)) * 0.18)
  )
}

cat("Saved rank5 downstream outputs\n")
