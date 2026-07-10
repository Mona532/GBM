library(SpaLinker)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_visium_tls_spalinker_official.R <results_dir>")
}

results_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
sample_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)

for (sample_dir in sample_dirs) {
  csv_path <- file.path(sample_dir, "tls_spot_scores_official_relaxed.csv")
  if (!file.exists(csv_path)) {
    next
  }
  message("[plot] ", basename(sample_dir))
  df <- read.csv(csv_path, check.names = FALSE)
  pos <- df[, c("x", "y")]

  SpotVisualize(
    pos = pos,
    meta = df[["TLS.score"]],
    title = paste0(basename(sample_dir), " TLS.score"),
    legend.name = "TLS.score",
    savefile = file.path(sample_dir, "tls_score_official.pdf"),
    scale_y_reverse = FALSE,
    p.width = 8,
    p.height = 7
  )

  SpotVisualize(
    pos = pos,
    meta = as.character(df[["TLS.region"]]),
    title = paste0(basename(sample_dir), " TLS.region"),
    legend.name = "TLS.region",
    cha.col = c(TLS = "#c32f27", nonTLS = "#ced4da"),
    savefile = file.path(sample_dir, "tls_region_official.pdf"),
    scale_y_reverse = FALSE,
    p.width = 8,
    p.height = 7
  )
}
