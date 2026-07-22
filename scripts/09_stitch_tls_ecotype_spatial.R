args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1) args[[1]] else "E:/GBM/results"
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

ecotypes <- sort(unique(weights$dominant_ecotype))
sig_map <- split(top_ct[top_ct$rank <= 3, "cell_type"], top_ct[top_ct$rank <= 3, "ecotype"])
score_csvs <- list.files(tls_root, pattern = "tls_spot_scores_official_relaxed.csv$", recursive = TRUE, full.names = TRUE)
score_values <- unlist(lapply(score_csvs, function(path) {
  df <- read.csv(path, check.names = FALSE)
  suppressWarnings(as.numeric(df$TLS.score))
}))
score_limits <- range(score_values, na.rm = TRUE, finite = TRUE)

read_spot_table <- function(sample_tag) {
  path <- file.path(tls_root, sample_tag, "tls_spot_scores_official_relaxed.csv")
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read.csv(path, check.names = FALSE)
  pos_cols <- if (all(c("x_pixel", "y_pixel") %in% colnames(df))) c("x_pixel", "y_pixel") else c("x", "y")
  names(df)[match(pos_cols, names(df))] <- c("x", "y")
  df
}

read_c2l_table <- function(sample_tag) {
  path <- file.path(c2l_root, sample_tag, "cell2loc_mean.csv")
  if (!file.exists(path)) {
    return(NULL)
  }
  c2l <- read.csv(path, row.names = 1, check.names = FALSE)
  c2l$barcode <- rownames(c2l)
  c2l
}

build_hulls <- function(df) {
  split_comp <- split(df, df$unit_id)
  hulls <- lapply(split_comp, function(sub) {
    if (nrow(sub) < 3) return(NULL)
    idx <- grDevices::chull(sub$x, sub$y)
    sub[idx, c("x", "y", "unit_id"), drop = FALSE]
  })
  hulls <- Filter(Negate(is.null), hulls)
  if (length(hulls) == 0) {
    return(NULL)
  }
  do.call(rbind, hulls)
}

build_small_components <- function(df) {
  small <- lapply(split(df, df$unit_id), function(sub) {
    if (nrow(sub) >= 3) return(NULL)
    sub[, c("x", "y", "unit_id"), drop = FALSE]
  })
  small <- Filter(Negate(is.null), small)
  if (length(small) == 0) {
    return(NULL)
  }
  do.call(rbind, small)
}

make_ecotype_plot <- function(sample_tag, ecotype, fill_limits = NULL) {
  spot_df <- read_spot_table(sample_tag)
  if (is.null(spot_df)) {
    return(NULL)
  }

  units <- weights$unit_id[weights$sample == sample_tag & weights$dominant_ecotype == ecotype]
  if (length(units) == 0) {
    return(NULL)
  }

  comp_spots <- spot_map[spot_map$unit_id %in% units, c("unit_id", "barcode"), drop = FALSE]
  if (nrow(comp_spots) == 0) {
    return(NULL)
  }

  target <- merge(spot_df, comp_spots, by = "barcode")
  if (nrow(target) == 0) {
    return(NULL)
  }

  c2l <- read_c2l_table(sample_tag)
  if (is.null(c2l)) {
    return(NULL)
  }
  sig_cells <- sig_map[[ecotype]]
  sig_cells <- sig_cells[sig_cells %in% colnames(c2l)]
  if (length(sig_cells) == 0) {
    return(NULL)
  }

  c2l$signature_abundance <- rowMeans(c2l[, sig_cells, drop = FALSE], na.rm = TRUE)
  target <- merge(target, c2l[, c("barcode", "signature_abundance"), drop = FALSE], by = "barcode", all.x = TRUE)
  target$tls_region_flag <- target$TLS.region == "TLS"

  bg <- spot_df[spot_df$in_tissue == 1, , drop = FALSE]
  hulls <- build_hulls(target)
  small <- build_small_components(target)
  sig_label <- paste(sig_map[[ecotype]], collapse = " + ")

  p <- ggplot() +
    geom_point(
      data = bg,
      aes(x = x, y = y),
      color = "grey88",
      size = 0.3
    ) +
    geom_point(
      data = target,
      aes(x = x, y = y, fill = signature_abundance, alpha = TLS.score),
      shape = 21,
      color = NA,
      size = 0.9,
      stroke = 0
    ) +
    geom_point(
      data = target[target$tls_region_flag, , drop = FALSE],
      aes(x = x, y = y),
      shape = 21,
      fill = NA,
      color = "#7f0000",
      stroke = 0.25,
      size = 1.0
    ) +
    coord_fixed() +
    scale_y_reverse() +
    scale_fill_gradientn(
      colors = c("#f7fcf5", "#a1d99b", "#238b45"),
      limits = fill_limits,
      name = "Top-cell abundance"
    ) +
    scale_alpha_continuous(
      range = c(0.20, 1),
      limits = score_limits,
      name = "TLS.score"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.subtitle = element_text(hjust = 0.5, size = 6),
      legend.position = "right"
    ) +
    ggtitle(sample_tag, subtitle = sig_label)

  if (!is.null(hulls) && nrow(hulls) > 0) {
    p <- p + geom_polygon(
      data = hulls,
      aes(x = x, y = y, group = unit_id),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.25
    )
  }
  if (!is.null(small) && nrow(small) > 0) {
    p <- p + geom_point(
      data = small,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      color = "black",
      size = 1.1,
      stroke = 0.3
    )
  }

  p
}

save_ecotype_batches <- function(ecotype, ncol = 5, nrow = 5) {
  sample_tags <- unique(weights$sample[weights$dominant_ecotype == ecotype])
  sample_tags <- sort(sample_tags)
  sig_cells <- sig_map[[ecotype]]
  fill_values <- unlist(lapply(sample_tags, function(sample_tag) {
    c2l <- read_c2l_table(sample_tag)
    if (is.null(c2l)) {
      return(numeric(0))
    }
    keep_cells <- sig_cells[sig_cells %in% colnames(c2l)]
    if (length(keep_cells) == 0) {
      return(numeric(0))
    }
    rowMeans(c2l[, keep_cells, drop = FALSE], na.rm = TRUE)
  }))
  fill_limits <- range(fill_values, na.rm = TRUE, finite = TRUE)
  plots <- lapply(sample_tags, make_ecotype_plot, ecotype = ecotype, fill_limits = fill_limits)
  keep <- vapply(plots, Negate(is.null), logical(1))
  plots <- plots[keep]
  sample_tags <- sample_tags[keep]
  if (length(plots) == 0) {
    return(invisible(NULL))
  }

  per_page <- ncol * nrow
  for (start in seq(1, length(plots), by = per_page)) {
    end <- min(start + per_page - 1, length(plots))
    idx <- start:end
    page <- wrap_plots(plots[idx], ncol = ncol, nrow = nrow) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )
    out <- file.path(
      tls_root,
      sprintf("tls_%s_spatial_grid_%02d-%02d.png", ecotype, start, end)
    )
    ggsave(out, page, width = 22, height = 20, dpi = 200, bg = "white")
    cat(sprintf("[%s] %s\n", ecotype, out))
  }
}

for (ecotype in ecotypes) {
  save_ecotype_batches(ecotype)
}
