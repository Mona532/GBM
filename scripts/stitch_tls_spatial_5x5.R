args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "E:/GBM/results/tls_core"
results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

sample_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
csvs <- file.path(sample_dirs, "tls_spot_scores_official_relaxed.csv")
csvs <- csvs[file.exists(csvs)]
tags <- basename(dirname(csvs))
spot_map <- read.csv(file.path(dirname(results_dir), "tls_component_spot_map.csv"), check.names = FALSE)
score_limits <- range(unlist(lapply(csvs, function(path) {
  df <- read.csv(path, check.names = FALSE)
  suppressWarnings(as.numeric(df$TLS.score))
})), na.rm = TRUE, finite = TRUE)

build_component_layers <- function(df, sample_tag) {
  comp_spots <- spot_map[spot_map$sample == sample_tag, c("barcode", "component_id"), drop = FALSE]
  if (nrow(comp_spots) == 0) {
    return(list(hulls = NULL, small = NULL))
  }
  merged <- merge(df[, c("barcode", "x", "y", "array_col", "array_row")], comp_spots, by = "barcode")
  if (nrow(merged) == 0) {
    return(list(hulls = NULL, small = NULL))
  }

  split_comp <- split(merged, merged$component_id)
  hulls <- lapply(split_comp, function(sub) {
    if (nrow(sub) < 3) return(NULL)
    idx <- grDevices::chull(sub$x, sub$y)
    sub[idx, c("x", "y", "component_id"), drop = FALSE]
  })
  hulls <- Filter(Negate(is.null), hulls)
  hulls <- if (length(hulls) > 0) do.call(rbind, hulls) else NULL

  small <- do.call(rbind, lapply(split_comp, function(sub) {
    if (nrow(sub) >= 3) return(NULL)
    sub[, c("x", "y", "component_id"), drop = FALSE]
  }))

  list(hulls = hulls, small = small)
}

make_tls_plot <- function(csv_path, mode = c("score", "region")) {
  mode <- match.arg(mode)
  df <- read.csv(csv_path, check.names = FALSE)
  pos_cols <- if (all(c("x_pixel", "y_pixel") %in% colnames(df))) c("x_pixel", "y_pixel") else c("x", "y")
  names(df)[match(pos_cols, names(df))] <- c("x", "y")
  sample_tag <- df$sample_id[[1]]
  comp_layers <- build_component_layers(df, sample_tag)

  p <- ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )

  if (mode == "score") {
    p <- p + geom_point(aes(color = TLS.score), size = 0.45) +
      scale_y_reverse() +
      scale_color_gradientn(
        colors = c("#0077b6", "lightyellow", "#c32f27"),
        limits = score_limits,
        name = "TLS.score"
      )
  } else {
    p <- p + geom_point(aes(color = TLS.region), size = 0.45) +
      scale_y_reverse() +
      scale_color_manual(
        values = c(TLS = "#c32f27", nonTLS = "#ced4da"),
        breaks = c("TLS", "nonTLS"),
        drop = FALSE,
        name = "TLS.region"
      )
  }

  if (!is.null(comp_layers$hulls) && nrow(comp_layers$hulls) > 0) {
    p <- p + geom_polygon(
      data = comp_layers$hulls,
      aes(x = x, y = y, group = component_id),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.25
    )
  }
  if (!is.null(comp_layers$small) && nrow(comp_layers$small) > 0) {
    p <- p + geom_point(
      data = comp_layers$small,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      stroke = 0.25,
      size = 0.8,
      fill = NA,
      color = "black"
    )
  }

  p + ggtitle(sample_tag)
}

save_batches <- function(mode = c("score", "region"), ncol = 5, nrow = 5) {
  mode <- match.arg(mode)
  per_page <- ncol * nrow
  for (start in seq(1, length(csvs), by = per_page)) {
    end <- min(start + per_page - 1, length(csvs))
    idx <- start:end
    plots <- lapply(csvs[idx], make_tls_plot, mode = mode)
    page <- wrap_plots(plots, ncol = ncol, nrow = nrow, guides = "collect") &
      theme(
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )
    out <- file.path(
      results_dir,
      sprintf("tls_%s_grid_%02d-%02d.png", mode, start, end)
    )
    ggsave(out, page, width = 22, height = 20, dpi = 200, bg = "white")
    cat(sprintf("[%s] %s\n", mode, out))
  }
}

save_batches("score")
save_batches("region")
