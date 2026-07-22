suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(ggplot2)
  library(msigdbr)
})

root <- "E:/GBM"
results_dir <- file.path(root, "results")
h5_dir <- file.path(root, "spatial_data_visium", "spatial_data_visium", "anndata_consolidated")
c2l_dir <- file.path(results_dir, "c2l_core_v1")

target_ecotype <- "E4"
ilc_cols <- c("ILC1", "ILC2", "ILC3")
all_axes <- c(ilc_cols, "ILC_total")
min_detect_rate <- 0.10
label_top_n <- 12L
bad_label_pattern <- "^(RPL|RPS|MT-|MTRNR|LINC|MALAT1$|NEAT1$)"
ora_min_genes <- 15L
ora_max_genes <- 500L
ora_top_terms <- 15L
neuro_pattern <- "(neuro|neuron|synap|axon|dendrit|glutam|gaba|seroton|dopamin|cholin|nervous|transmission|vesicle)"

read_h5_counts <- function(path) {
  genes <- h5read(path, "var/_index")
  barcodes <- h5read(path, "obs/_index")
  xdata <- as.numeric(h5read(path, "X/data"))
  idx <- as.integer(h5read(path, "X/indices")) + 1L
  indptr <- as.integer(h5read(path, "X/indptr"))
  row_id <- rep(seq_len(length(barcodes)), diff(indptr))

  mat_obs_gene <- sparseMatrix(
    i = row_id,
    j = idx,
    x = xdata,
    dims = c(length(barcodes), length(genes))
  )
  rownames(mat_obs_gene) <- barcodes
  colnames(mat_obs_gene) <- genes

  feature_cats <- h5read(path, "var/feature_types/categories")
  feature_codes <- as.integer(h5read(path, "var/feature_types/codes"))
  feature_types <- feature_cats[feature_codes + 1L]

  list(
    mat_counts = mat_obs_gene,
    barcodes = barcodes,
    genes = genes,
    feature_types = feature_types
  )
}

log1p_cpm <- function(mat) {
  lib <- Matrix::rowSums(mat)
  lib[lib == 0] <- 1
  scale_vec <- 10000 / lib
  norm <- Diagonal(x = scale_vec) %*% mat
  log1p(as.matrix(norm))
}

make_corr_table <- function(expr_mat, axis_vals, axis_name) {
  rho <- suppressWarnings(cor(t(expr_mat), axis_vals, method = "spearman"))
  rho <- as.numeric(rho)
  n <- ncol(expr_mat)
  rho_clip <- pmin(pmax(rho, -0.999999), 0.999999)
  t_stat <- rho_clip * sqrt((n - 2) / pmax(1e-12, 1 - rho_clip^2))
  pval <- 2 * pt(-abs(t_stat), df = n - 2)
  out <- data.frame(
    gene = rownames(expr_mat),
    rho = rho,
    pvalue = pval,
    mean_expr = rowMeans(expr_mat),
    detect_rate = rowMeans(expr_mat > 0),
    axis = axis_name,
    stringsAsFactors = FALSE
  )
  out$fdr <- p.adjust(out$pvalue, method = "BH")
  out <- out[order(-out$rho, out$pvalue), ]
  rownames(out) <- NULL
  out
}

run_simple_ora <- function(genes_of_interest, background_genes, term2gene, direction, axis_name) {
  genes_of_interest <- unique(intersect(genes_of_interest, background_genes))
  background_genes <- unique(background_genes)
  if (length(genes_of_interest) < ora_min_genes) {
    return(data.frame())
  }

  term_split <- split(term2gene$gene_symbol, paste(term2gene$gs_collection, term2gene$gs_name, sep = "::"))
  n <- length(genes_of_interest)
  N <- length(background_genes)
  out <- lapply(names(term_split), function(term_id) {
    term_genes <- unique(intersect(term_split[[term_id]], background_genes))
    k <- length(term_genes)
    x <- length(intersect(genes_of_interest, term_genes))
    if (k < 5 || x < 2) return(NULL)
    p <- phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE)
    pieces <- strsplit(term_id, "::", fixed = TRUE)[[1]]
    data.frame(
      axis = axis_name,
      direction = direction,
      gene_set = pieces[1],
      term = pieces[2],
      overlap = x,
      term_size = k,
      input_size = n,
      background_size = N,
      gene_ratio = x / n,
      pvalue = p,
      hit_genes = paste(sort(intersect(genes_of_interest, term_genes)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out$fdr <- p.adjust(out$pvalue, method = "BH")
  out$neuro_related <- grepl(neuro_pattern, out$term, ignore.case = TRUE)
  out <- out[order(out$fdr, -out$overlap, out$pvalue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_ora_summary <- function(df, axis_name, out_jpg, out_pdf) {
  if (nrow(df) == 0) return(invisible(NULL))
  keep <- do.call(rbind, lapply(split(df, interaction(df$direction, df$gene_set, drop = TRUE)), function(x) head(x, ora_top_terms)))
  keep$term <- factor(keep$term, levels = rev(unique(keep$term)))
  keep$score <- -log10(pmax(keep$fdr, 1e-300))

  p <- ggplot(keep, aes(x = score, y = term, size = overlap, color = neuro_related)) +
    geom_point(alpha = 0.85) +
    facet_grid(direction ~ gene_set, scales = "free_y", space = "free_y") +
    scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
    labs(
      title = paste0(axis_name, " pathway enrichment"),
      x = "-log10(FDR)",
      y = NULL,
      color = "Neuro-related"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(out_jpg, p, width = 10.5, height = 7.5, dpi = 300)
  ggsave(out_pdf, p, width = 10.5, height = 7.5)
}

plot_corr_summary <- function(df, axis_name, out_jpg, out_pdf) {
  df$neglog10p <- -log10(pmax(df$pvalue, 1e-300))
  df$sig <- df$fdr < 0.05

  label_df <- df[!grepl(bad_label_pattern, df$gene, ignore.case = TRUE), , drop = FALSE]
  top_pos <- head(label_df[order(-label_df$rho, label_df$pvalue), ], label_top_n)
  top_neg <- head(label_df[order(label_df$rho, label_df$pvalue), ], label_top_n)
  lab <- unique(rbind(top_pos, top_neg))

  p <- ggplot(df, aes(x = rho, y = neglog10p)) +
    geom_point(data = df[!df$sig, , drop = FALSE], color = "#C7CCD6", size = 1.1, alpha = 0.7) +
    geom_point(data = df[df$sig, , drop = FALSE], color = "#B2182B", size = 1.2, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "#666666") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "#666666") +
    geom_text(data = lab, aes(label = gene), size = 2.3, vjust = -0.2, check_overlap = TRUE) +
    labs(
      title = paste0("E4 spot-level gene correlation with ", axis_name),
      x = paste0("Spearman rho with ", axis_name, " abundance"),
      y = "-log10(p-value)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(out_jpg, p, width = 7.2, height = 5.8, dpi = 300)
  ggsave(out_pdf, p, width = 7.2, height = 5.8)
}

weights <- read.csv(file.path(results_dir, "tls_compnmf_rank5_unit_weights.csv"), check.names = FALSE, stringsAsFactors = FALSE)
spot_map <- read.csv(file.path(results_dir, "tls_component_spot_map.csv"), check.names = FALSE, stringsAsFactors = FALSE)
target_units <- weights$unit_id[weights$dominant_ecotype == target_ecotype]
e4_spots <- spot_map[spot_map$unit_id %in% target_units, c("sample", "unit_id", "barcode"), drop = FALSE]
if (nrow(e4_spots) == 0) stop("No E4 component spots found.")

expr_list <- list()
abund_list <- list()
meta_list <- list()

for (sid in sort(unique(e4_spots$sample))) {
  h5_path <- file.path(h5_dir, paste0(sid, ".h5ad"))
  c2l_csv <- file.path(c2l_dir, sid, "cell2loc_mean.csv")
  if (!file.exists(h5_path) || !file.exists(c2l_csv)) next

  sample_spots <- e4_spots[e4_spots$sample == sid, , drop = FALSE]
  target_barcodes <- unique(sample_spots$barcode)

  dat <- read_h5_counts(h5_path)
  gene_mask <- dat$feature_types == "Gene Expression"
  gene_idx <- which(gene_mask)
  mat_counts <- dat[["mat_counts"]]
  mat_counts <- mat_counts[, gene_idx, drop = FALSE]

  c2l <- read.csv(c2l_csv, row.names = 1, check.names = FALSE)
  missing_expr <- setdiff(target_barcodes, rownames(mat_counts))
  missing_c2l <- setdiff(target_barcodes, rownames(c2l))
  if (length(missing_expr) > 0 || length(missing_c2l) > 0) {
    msg <- c(
      sprintf("Barcode validation failed for sample %s", sid),
      if (length(missing_expr) > 0) {
        sprintf(
          "Missing in h5ad (%d): %s",
          length(missing_expr),
          paste(head(missing_expr, 10), collapse = ", ")
        )
      },
      if (length(missing_c2l) > 0) {
        sprintf(
          "Missing in cell2loc (%d): %s",
          length(missing_c2l),
          paste(head(missing_c2l, 10), collapse = ", ")
        )
      }
    )
    stop(paste(msg, collapse = "\n"))
  }

  mat_counts <- mat_counts[target_barcodes, , drop = FALSE]
  expr <- log1p_cpm(mat_counts)
  c2l <- c2l[target_barcodes, ilc_cols, drop = FALSE]

  expr_list[[sid]] <- expr
  abund_list[[sid]] <- as.matrix(c2l)
  meta_list[[sid]] <- sample_spots[match(target_barcodes, sample_spots$barcode), , drop = FALSE]
}

if (length(expr_list) == 0) stop("No E4 sample expression loaded.")

common_genes <- Reduce(intersect, lapply(expr_list, colnames))
expr_all <- do.call(rbind, lapply(expr_list, function(x) x[, common_genes, drop = FALSE]))
abund_all <- do.call(rbind, abund_list)
meta_all <- do.call(rbind, meta_list)
rownames(meta_all) <- meta_all$barcode
meta_all <- meta_all[rownames(expr_all), , drop = FALSE]
abund_all <- as.data.frame(abund_all)
abund_all$ILC_total <- rowSums(abund_all[, ilc_cols, drop = FALSE])

detect_rate <- colMeans(expr_all > 0)
keep_genes <- names(detect_rate)[detect_rate >= min_detect_rate]
expr_all <- expr_all[, keep_genes, drop = FALSE]
expr_gene_by_spot <- t(expr_all)

summary_df <- data.frame(
  ecotype = target_ecotype,
  n_units = length(unique(e4_spots$unit_id)),
  n_samples = length(unique(e4_spots$sample)),
  n_spots = nrow(expr_all),
  n_genes_tested = nrow(expr_gene_by_spot),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(results_dir, "e4_spot_ilc_gene_corr_summary.csv"), row.names = FALSE)
write.csv(meta_all, file.path(results_dir, "e4_spot_ilc_gene_corr_spot_metadata.csv"), row.names = FALSE)

for (axis_name in all_axes) {
  axis_vals <- abund_all[, axis_name]
  corr_df <- make_corr_table(expr_gene_by_spot, axis_vals, axis_name)
  write.csv(corr_df, file.path(results_dir, paste0("e4_spot_gene_corr_", axis_name, ".csv")), row.names = FALSE)

  label_df <- corr_df[!grepl(bad_label_pattern, corr_df$gene, ignore.case = TRUE), , drop = FALSE]
  top_tbl <- rbind(
    cbind(direction = "positive", head(label_df[order(-label_df$rho, label_df$pvalue), ], 30)),
    cbind(direction = "negative", head(label_df[order(label_df$rho, label_df$pvalue), ], 30))
  )
  write.csv(top_tbl, file.path(results_dir, paste0("e4_spot_gene_corr_", axis_name, "_top30_posneg.csv")), row.names = FALSE)

  plot_corr_summary(
    corr_df,
    axis_name,
    file.path(results_dir, paste0("fig_e4_spot_gene_corr_", axis_name, ".jpg")),
    file.path(results_dir, paste0("fig_e4_spot_gene_corr_", axis_name, ".pdf"))
  )

  cat(sprintf("%s: wrote %d genes\n", axis_name, nrow(corr_df)))
}

term2gene <- msigdbr(species = "Homo sapiens") |>
  subset((gs_collection == "C2" & gs_subcollection == "CP:REACTOME") | gs_collection == "C5") |>
  subset(gs_collection != "C5" | gs_subcollection == "GO:BP") |>
  unique()

ilc_total_df <- read.csv(file.path(results_dir, "e4_spot_gene_corr_ILC_total.csv"), check.names = FALSE, stringsAsFactors = FALSE)
bg_genes <- ilc_total_df$gene
pos_genes <- subset(
  ilc_total_df,
  rho > 0 & fdr < 0.05 & !grepl(bad_label_pattern, gene, ignore.case = TRUE)
)$gene
neg_genes <- subset(
  ilc_total_df,
  rho < 0 & fdr < 0.05 & !grepl(bad_label_pattern, gene, ignore.case = TRUE)
)$gene

if (length(pos_genes) > ora_max_genes) {
  pos_ranked <- subset(ilc_total_df, gene %in% pos_genes)
  pos_genes <- head(pos_ranked[order(-pos_ranked$rho, pos_ranked$pvalue), "gene"], ora_max_genes)
}
if (length(neg_genes) > ora_max_genes) {
  neg_ranked <- subset(ilc_total_df, gene %in% neg_genes)
  neg_genes <- head(neg_ranked[order(neg_ranked$rho, neg_ranked$pvalue), "gene"], ora_max_genes)
}

ora_df <- rbind(
  run_simple_ora(pos_genes, bg_genes, term2gene, "positive", "ILC_total"),
  run_simple_ora(neg_genes, bg_genes, term2gene, "negative", "ILC_total")
)

if (nrow(ora_df) > 0) {
  write.csv(ora_df, file.path(results_dir, "e4_spot_gene_corr_ILC_total_pathway_enrichment.csv"), row.names = FALSE)
  write.csv(head(subset(ora_df, neuro_related), 50), file.path(results_dir, "e4_spot_gene_corr_ILC_total_pathway_neuro_hits_top50.csv"), row.names = FALSE)
  plot_ora_summary(
    ora_df,
    "ILC_total",
    file.path(results_dir, "fig_e4_spot_gene_corr_ILC_total_pathway.jpg"),
    file.path(results_dir, "fig_e4_spot_gene_corr_ILC_total_pathway.pdf")
  )
  cat(sprintf("ILC_total ORA: wrote %d terms\n", nrow(ora_df)))
} else {
  cat("ILC_total ORA: no enriched terms written\n")
}

cat("Saved E4 spot-level ILC gene correlation outputs\n")
print(summary_df)
