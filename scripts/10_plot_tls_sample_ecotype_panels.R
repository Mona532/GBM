args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1) args[[1]] else "E:/GBM/results"
top_n <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
results_root <- normalizePath(results_root, winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

weights <- read.csv(file.path(results_root, "tls_compnmf_rank5_unit_weights.csv"), check.names = FALSE)
spot_map <- read.csv(file.path(results_root, "tls_component_spot_map.csv"), check.names = FALSE)
top_ct <- read.csv(file.path(results_root, "tls_compnmf_rank5_top_celltypes.csv"), check.names = FALSE)

tls_root <- file.path(results_root, "tls_core")
c2l_root <- file.path(results_root, "c2l_core_v1")
out_dir <- file.path(tls_root, "sample_ecotype_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ecotypes <- sort(unique(weights$dominant_ecotype))
score_csvs <- list.files(tls_root, pattern = "tls_spot_scores_official_relaxed.csv$", recursive = TRUE, full.names = TRUE)
score_values <- unlist(lapply(score_csvs, function(path) {
  df <- read.csv(path, check.names = FALSE)
  suppressWarnings(as.numeric(df$TLS.score))
}))
score_limits <- range(score_values, na.rm = TRUE, finite = TRUE)

read_spot_table <- function(sample_tag) {
  path <- file.path(tls_root, sample_tag, "tls_spot_scores_official_relaxed.csv")
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, check.names = FALSE)
  pos_cols <- if (all(c("x_pixel", "y_pixel") %in% colnames(df))) c("x_pixel", "y_pixel") else c("x", "y")
  names(df)[match(pos_cols, names(df))] <- c("x", "y")
  df
}

read_c2l_table <- function(sample_tag) {
  path <- file.path(c2l_root, sample_tag, "cell2loc_mean.csv")
  if (!file.exists(path)) return(NULL)
  c2l <- read.csv(path, row.names = 1, check.names = FALSE)
  c2l$barcode <- rownames(c2l)
  c2l
}

panel_base <- function(full_df) {
  ggplot(full_df, aes(x = x, y = y)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      legend.position = "right"
    )
}

make_score_panel <- function(full_df, title_text) {
  panel_base(full_df) +
    geom_point(aes(color = TLS.score), size = 1.2) +
    scale_color_gradientn(
      colors = c("#0077b6", "lightyellow", "#c32f27"),
      limits = score_limits,
      name = "TLS.score"
    ) +
    ggtitle(title_text)
}

make_region_panel <- function(full_df, target_df) {
  panel_base(full_df) +
    geom_point(color = "#ced4da", size = 1.0) +
    geom_point(
      data = target_df,
      aes(x = x, y = y, color = TLS.region),
      size = 1.2,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(TLS = "#c32f27", nonTLS = "#ced4da"),
      breaks = c("TLS", "nonTLS"),
      drop = FALSE,
      name = "TLS.region"
    ) +
    ggtitle("TLS.region")
}

make_abundance_panel <- function(full_df, target_df, cell_type) {
  vals <- target_df[[cell_type]]
  lims <- range(vals, na.rm = TRUE, finite = TRUE)
  if (!all(is.finite(lims)) || diff(lims) == 0) {
    lims <- c(min(vals, na.rm = TRUE), min(vals, na.rm = TRUE) + 1e-6)
  }
  panel_base(full_df) +
    geom_point(color = "grey88", size = 1.0) +
    geom_point(
      data = target_df,
      aes(x = x, y = y, color = .data[[cell_type]]),
      size = 1.2,
      inherit.aes = FALSE
    ) +
    scale_color_viridis_c(
      option = "magma",
      limits = lims,
      name = cell_type
    ) +
    ggtitle(cell_type)
}

sample_ecotype_summary <- aggregate(n_spots ~ sample + dominant_ecotype, data = weights, FUN = sum)

for (ecotype in ecotypes) {
  sample_sub <- sample_ecotype_summary[sample_ecotype_summary$dominant_ecotype == ecotype, , drop = FALSE]
  sample_sub <- sample_sub[order(-sample_sub$n_spots, sample_sub$sample), , drop = FALSE]
  sample_tags <- head(sample_sub$sample, top_n)
  cell_types <- top_ct$cell_type[top_ct$ecotype == ecotype & top_ct$rank <= 5]

  for (sample_tag in sample_tags) {
    spot_df <- read_spot_table(sample_tag)
    c2l <- read_c2l_table(sample_tag)
    weight_sub <- weights[weights$sample == sample_tag & weights$dominant_ecotype == ecotype, , drop = FALSE]
    units <- weight_sub$unit_id
    if (is.null(spot_df) || is.null(c2l) || length(units) == 0) next

    top_unit <- weight_sub$unit_id[which.max(weight_sub$n_spots)]
    comp_spots <- spot_map[spot_map$unit_id == top_unit, c("unit_id", "barcode"), drop = FALSE]
    if (nrow(comp_spots) == 0) next

    target <- merge(spot_df, comp_spots, by = "barcode")
    if (nrow(target) == 0) next

    keep_cells <- cell_types[cell_types %in% colnames(c2l)]
    if (length(keep_cells) == 0) next
    target <- merge(target, c2l[, c("barcode", keep_cells), drop = FALSE], by = "barcode", all.x = TRUE)
    full_df <- merge(
      spot_df[spot_df$in_tissue == 1, , drop = FALSE],
      c2l[, c("barcode", keep_cells), drop = FALSE],
      by = "barcode",
      all.x = TRUE
    )

    panels <- c(
      list(make_score_panel(full_df, paste0(ecotype, " TLS.score"))),
      list(make_region_panel(full_df, target)),
      lapply(keep_cells, function(ct) make_abundance_panel(full_df, target, ct))
    )

    page <- wrap_plots(panels, ncol = 3) &
      theme(
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    out <- file.path(out_dir, sprintf("%s__%s__%s_panels.png", ecotype, sample_tag, gsub("__", "_", top_unit, fixed = TRUE)))
    ggsave(out, page, width = 16, height = 10, dpi = 220, bg = "white")
    cat(sprintf("[%s] %s\n", ecotype, out))
  }
}
