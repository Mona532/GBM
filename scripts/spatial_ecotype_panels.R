# ============================================================
# spatial_ecotype_panels.R — 5 figures, one per ecotype,
# each sample one row: TLS score → ecotype → celltype1 → celltype2
# ============================================================
library(rhdf5); library(SpaLinker); library(cowplot)

root <- "E:/GBM"; out_dir <- file.path(root, "results", "spatial_ecotype_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spot_map <- read.csv(file.path(root, "results", "tls_component_spot_map.csv"))
weights  <- read.csv(file.path(root, "results", "tls_compnmf_rank5_unit_weights.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
eco_label <- setNames(weights$dominant_ecotype, weights$unit_id)
spot_map$ecotype <- eco_label[spot_map$unit_id]

# defining cell types per ecotype
eco_cts <- list(
  E1 = c("B",        "Tfh-like_CD4"),
  E2 = c("HEV-like_endothelial", "Glioma"),
  E3 = c("Monocyte", "Macrophage"),
  E4 = c("ILC3",     "ILC2"),
  E5 = c("Glial",    "Tfh-like_CD4")
)

read_h5 <- function(path) {
  obs_names <- as.character(h5read(path, "/obs/_index"))
  coords <- t(h5read(path, "/obsm/spatial"))
  colnames(coords) <- c("x", "y"); rownames(coords) <- obs_names
  q05 <- t(h5read(path, "/obsm/c2l_ilc_q05"))
  ct   <- h5read(path, "/uns/c2l_ilc_cell_types")
  colnames(q05) <- ct; rownames(q05) <- obs_names
  list(coords = as.data.frame(coords), q05 = q05, ct = ct)
}

for (eco in paste0("E", 1:5)) {
  eco_spots <- spot_map[spot_map$ecotype == eco & !is.na(spot_map$ecotype), ]
  samples   <- sort(unique(eco_spots$sample))
  cts <- eco_cts[[eco]]
  plots <- list()

  for (sid in samples) {
    h5 <- file.path(root, "spatial_data_visium", "spatial_data_visium",
                    "anndata_core", paste0(sid, ".h5ad"))
    if (!file.exists(h5)) next

    dat <- read_h5(h5)
    ss  <- eco_spots[eco_spots$sample == sid, ]
    n_spots <- nrow(ss)
    n_comp  <- length(unique(ss$component_id))

    # panel 1: TLS score
    tls <- tryCatch(
      read.csv(file.path(root, "results", "tls_core", sid,
               "tls_spot_scores_official_relaxed.csv"), row.names = 1),
      error = function(e) NULL)
    tls_score <- rep(0, nrow(dat$coords))
    if (!is.null(tls)) {
      tls$TLS.score[is.na(tls$TLS.score)] <- 0
      tls_score <- tls$TLS.score
    }
    p1 <- SpotVisualize(pos = dat$coords, meta = tls_score, title = "TLS",
          size = 1.2, f.color = c("#0077b6", "lightyellow", "#c32f27"),
          p.width = 1.8, p.height = 1.8, return = TRUE)

    eco_pos <- rep(NA_character_, nrow(dat$coords))
    eco_pos[match(ss$barcode, rownames(dat$coords))] <- eco
    p2 <- SpotVisualize(pos = dat$coords, meta = eco_pos,
          title = eco, size = 1.5, p.width = 1.8, p.height = 1.8,
          return = TRUE, na.col = "lightgrey")

    p_extra <- list()
    for (ct in cts) {
      if (!ct %in% dat$ct) next
      vals <- rep(NA_real_, nrow(dat$coords))
      vals[match(ss$barcode, rownames(dat$coords))] <- dat$q05[ss$barcode, ct]
      p_extra[[length(p_extra) + 1]] <- SpotVisualize(pos = dat$coords,
        meta = vals, title = ct, size = 1.2,
        f.color = c("#0077b6", "lightyellow", "#c32f27"),
        p.width = 1.8, p.height = 1.8, return = TRUE, na.col = "lightgrey")
    }

    row_plots <- c(list(p1, p2), p_extra)
    plots[[length(plots) + 1]] <- plot_grid(
      plotlist = row_plots, ncol = length(row_plots),
      labels = sid, label_size = 6, label_fontface = "plain",
      hjust = -0.05)
  }

  if (length(plots) == 0) next
  ncol <- 1
  combined <- plot_grid(plotlist = plots, ncol = ncol)
  w <- (2 + length(cts)) * 2.0; h <- length(plots) * 2.0
  pdf(file.path(out_dir, paste0("fig_spatial_", eco, "_panels.pdf")),
      width = w, height = h)
  print(combined); dev.off()
  cat(sprintf("%s: %d samples\n", eco, length(samples)))
}
cat("Done\n")
