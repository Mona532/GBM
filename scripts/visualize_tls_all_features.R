# Visualize all TLS features per sample using SpaLinker::SpotVisualize
library(SpaLinker)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "E:/GBM/results/tls_official_cut01"
results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)

sample_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
cat(sprintf("Found %d samples\n", length(sample_dirs)))

for (sample_dir in sample_dirs) {
  csv_path <- file.path(sample_dir, "tls_spot_scores_official_relaxed.csv")
  if (!file.exists(csv_path)) next

  tag <- basename(sample_dir)
  cat(sprintf("[viz] %s\n", tag))
  df <- read.csv(csv_path, check.names = FALSE)

  # Use pixel coords for visualization (looks like the tissue)
  pos <- df[, c("x_pixel", "y_pixel")]
  colnames(pos) <- c("x", "y")

  # Panel 1: TLS score
  SpotVisualize(
    pos = pos,
    meta = df$TLS.score,
    title = paste0(tag, " | TLS.score"),
    legend.name = "TLS.score",
    savefile = file.path(sample_dir, "tls_score_official.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 2: TLS region
  SpotVisualize(
    pos = pos,
    meta = as.character(df$TLS.region),
    title = paste0(tag, " | TLS.region"),
    legend.name = "TLS.region",
    cha.col = c(TLS = "#c32f27", nonTLS = "#ced4da"),
    savefile = file.path(sample_dir, "tls_region_official.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 3: Plasma/B.cells abundance
  SpotVisualize(
    pos = pos,
    meta = df$Plasma_B_cells,
    title = paste0(tag, " | Plasma/B.cells"),
    legend.name = "Abundance",
    f.color = c("#f7fbff", "#2171b5"),
    savefile = file.path(sample_dir, "plasma_b_cells.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 4: T.cells abundance
  SpotVisualize(
    pos = pos,
    meta = df$T_cells,
    title = paste0(tag, " | T.cells"),
    legend.name = "Abundance",
    f.color = c("#fff5f0", "#cb181d"),
    savefile = file.path(sample_dir, "t_cells.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 5: B-T co-distribution
  SpotVisualize(
    pos = pos,
    meta = df$Plasma_B.cells_T.cells,
    title = paste0(tag, " | B-T co-distribution"),
    legend.name = "Co-dist",
    f.color = c("#f7fcf5", "#238b45"),
    savefile = file.path(sample_dir, "bt_codistribution.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 6: LC.50sig enrichment
  SpotVisualize(
    pos = pos,
    meta = df$LC.50sig,
    title = paste0(tag, " | LC.50sig"),
    legend.name = "Score",
    f.color = c("#f7f4f9", "#6a51a3"),
    savefile = file.path(sample_dir, "lc50_sig.pdf"),
    p.width = 8, p.height = 7
  )

  # Panel 7: Dominant niche cluster
  niche_vals <- df$dominant_niche_cluster
  SpotVisualize(
    pos = pos,
    meta = niche_vals,
    title = paste0(tag, " | Spatial niche"),
    legend.name = "Niche",
    savefile = file.path(sample_dir, "spatial_niche.pdf"),
    p.width = 8, p.height = 7
  )
}

cat("Done\n")
