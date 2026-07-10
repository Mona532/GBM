# =============================================================================
# SpaLinker TLS Detection — CORRECT workflow per paper (Cell Genomics, 2025)
#
# Three TLS features → geometric mean:
#   1. LC.50sig enrichment (molecular TLS signature)
#   2. Plasma/B cell × T cell co-distribution (cellular co-localization)
#   3. (additional feature if available)
#
# Steps:
#   1. Load + preprocess Visium data
#   2. Define B/Plasma + T cell gene signatures → score each spot
#   3. Calculate B_Plasma-T cell co-distribution (product of normalized scores)
#   4. Calculate LC.50sig enrichment
#   5. Get spatial domains (BayesSpace via BayesCluster or Seurat clusters)
#   6. CalTLSfea(data=[LC.50sig, B_T_codist], st_pos, cluster=domains)
# =============================================================================

library(Seurat)
library(SpaLinker)

# ---- B/Plasma cell gene signature ----
b_plasma_genes <- c(
  "MS4A1", "CD19", "CD79A", "CD79B", "CD22", "CD37", "CD40",
  "PAX5", "BLK", "BANK1", "CD27", "TNFRSF13C", "CXCR5",
  "SDC1", "MZB1", "TNFRSF17", "XBP1", "IRF4", "JCHAIN",
  "DERL3", "SLAMF7", "IGHM", "IGHD", "IGHG1", "IGHA1",
  "IGKC", "IGLC2", "FCRL5", "POU2AF1", "SPIB"
)

# ---- T cell gene signature ----
t_cell_genes <- c(
  "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B",
  "TRAC", "TRBC1", "TRBC2", "CD28", "CD247", "LCK", "ZAP70",
  "ITK", "TXNIP", "CCL5", "GZMK", "GZMA", "GZMH", "NKG7"
)

# ---- LC.50sig (from CuratedSig.lt) ----
lc50_genes <- CuratedSig.lt$Immune$TLS$LC.50sig

# ---- Sample metadata (IDH-WT GBM) ----
samples <- list(
  list(gsm = "GSM7596587", name_tag = "mgh258",   region = "unannotated",  mgmt = "methylated"),
  list(gsm = "GSM7596588", name_tag = "zh881inf",  region = "infiltrating", mgmt = "partial"),
  list(gsm = "GSM7596589", name_tag = "zh881t1",   region = "T1",           mgmt = "partial"),
  list(gsm = "GSM7596590", name_tag = "zh916bulk", region = "bulk",         mgmt = "methylated"),
  list(gsm = "GSM7596591", name_tag = "zh916inf",  region = "infiltrating", mgmt = "methylated"),
  list(gsm = "GSM7596592", name_tag = "zh916t1",   region = "T1",           mgmt = "methylated"),
  list(gsm = "GSM7596593", name_tag = "zh1007inf", region = "infiltrating", mgmt = "NA"),
  list(gsm = "GSM7596594", name_tag = "zh1007nec", region = "necrotic",     mgmt = "NA"),
  list(gsm = "GSM7596595", name_tag = "zh1019inf", region = "infiltrating", mgmt = "NA"),
  list(gsm = "GSM7596596", name_tag = "zh1019t1",  region = "T1",           mgmt = "NA"),
  list(gsm = "GSM7596597", name_tag = "zh8811a",   region = "bulk",         mgmt = "partial"),
  list(gsm = "GSM7596598", name_tag = "zh8811b",   region = "bulk",         mgmt = "partial"),
  list(gsm = "GSM7596599", name_tag = "zh8812",    region = "bulk",         mgmt = "partial")
)
# Run all 13 GBM samples

DATA_DIR <- "E:/GBM/GSE237183_RAW"
OUT_DIR  <- "E:/GBM/results/spalinker_tls_v3"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results_list <- list()

for (i in seq_along(samples)) {
  samp <- samples[[i]]
  cat(sprintf("\n========== [%d/%d] %s ==========\n", i, length(samples), samp$name_tag))

  file_prefix <- file.path(DATA_DIR, paste0(samp$gsm, "_", samp$name_tag))
  h5_file <- paste0(file_prefix, "_filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_file)) { cat("  SKIPPING: no h5\n"); next }

  # ---- Setup temp Visium directory ----
  tmp_dir <- file.path(OUT_DIR, "tmp", samp$gsm)
  spatial_subdir <- file.path(tmp_dir, "spatial")
  dir.create(spatial_subdir, recursive = TRUE, showWarnings = FALSE)
  file.copy(h5_file, file.path(tmp_dir, "filtered_feature_bc_matrix.h5"), overwrite = TRUE)
  for (ext in c("tissue_positions_list.csv.gz", "scalefactors_json.json.gz",
                "tissue_lowres_image.png.gz")) {
    src <- paste0(file_prefix, "_", ext)
    dest_name <- sub("\\.gz$", "", ext)
    if (file.exists(src)) {
      R.utils::gunzip(src, destname = file.path(spatial_subdir, dest_name),
                      overwrite = TRUE, remove = FALSE)
    }
  }

  tryCatch({
    # ---- Step 1: Load + Preprocess ----
    cat("  [1] Loading Visium...\n")
    se <- Load10X_Spatial(data.dir = tmp_dir, assay = "Spatial",
                          filter.matrix = TRUE, slice = samp$name_tag)
    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)
    cat(sprintf("    Spots: %d\n", ncol(se)))
    if (ncol(se) < 50) next

    cat("  [2] SCTransform + clustering...\n")
    se <- SCTransform(se, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)
    se <- RunPCA(se, assay = "SCT", verbose = FALSE)
    se <- FindNeighbors(se, dims = 1:20, verbose = FALSE)
    se <- FindClusters(se, resolution = 0.6, verbose = FALSE)

    # ---- Step 2: Score B/Plasma cells and T cells per spot ----
    cat("  [3] Scoring B/Plasma + T cell signatures...\n")
    expr <- GetAssayData(se, assay = "Spatial", slot = "counts")

    b_avail <- intersect(b_plasma_genes, rownames(expr))
    t_avail <- intersect(t_cell_genes, rownames(expr))
    lc50_avail <- intersect(lc50_genes, rownames(expr))
    cat(sprintf("    B/Plasma genes: %d/%d, T genes: %d/%d, LC.50sig: %d/%d\n",
                length(b_avail), length(b_plasma_genes),
                length(t_avail), length(t_cell_genes),
                length(lc50_avail), length(lc50_genes)))

    cell_scores <- GsetScore(expr = expr,
      geneset = list(B_Plasma = b_avail, T_cell = t_avail),
      method  = "AddModuleScore", scale = FALSE, verbose = FALSE)

    # ---- Step 3: B/Plasma × T cell co-distribution ----
    cat("  [4] Computing B_Plasma × T cell co-distribution...\n")
    # Normalize each to [0,1]
    norm_bp  <- (cell_scores$B_Plasma - min(cell_scores$B_Plasma)) /
                (max(cell_scores$B_Plasma) - min(cell_scores$B_Plasma))
    norm_t   <- (cell_scores$T_cell - min(cell_scores$T_cell)) /
                (max(cell_scores$T_cell) - min(cell_scores$T_cell))
    b_t_codist <- norm_bp * norm_t
    cat(sprintf("    Co-dist range: [%.4f, %.4f]\n", min(b_t_codist), max(b_t_codist)))

    # ---- Step 4: LC.50sig enrichment ----
    cat("  [5] Computing LC.50sig enrichment...\n")
    lc50_score <- GsetScore(expr = expr,
      geneset = list(LC50sig = lc50_avail),
      method  = "AddModuleScore", scale = FALSE, verbose = FALSE)$LC50sig

    # ---- Step 5: Build feature matrix for CalTLSfea ----
    cat("  [6] Building TLS feature matrix...\n")
    tls_features <- data.frame(
      LC50sig     = lc50_score,
      B_T_codist  = b_t_codist,
      row.names   = colnames(se),
      check.names = FALSE
    )
    cat(sprintf("    Features: %d spots × %d features\n", nrow(tls_features), ncol(tls_features)))

    # Spatial positions
    st_pos <- se@images[[samp$name_tag]]@coordinates[, c("imagecol", "imagerow")]
    colnames(st_pos) <- c("x", "y")
    rownames(st_pos) <- colnames(se)

    # Domain labels (for neighborhood constraint)
    domains <- as.character(Idents(se))
    names(domains) <- colnames(se)

    # ---- Step 6: Run CalTLSfea ----
    cat("  [7] Running CalTLSfea...\n")
    tls_result <- CalTLSfea(
      data       = tls_features,
      st_pos     = st_pos,
      cluster    = domains,
      cutoff     = 0.2,
      filt.dist  = 4,
      filt.spots = 2
    )

    # ---- Results ----
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]
    se$B_Plasma   <- cell_scores$B_Plasma[colnames(se)]
    se$T_cell     <- cell_scores$T_cell[colnames(se)]
    se$BT_codist  <- b_t_codist[colnames(se)]
    se$LC50sig    <- lc50_score[colnames(se)]

    n_tls   <- sum(se$TLS_region == "TLS", na.rm = TRUE)
    n_total <- ncol(se)
    cat(sprintf("    TLS spots: %d / %d (%.1f%%)\n", n_tls, n_total, 100*n_tls/n_total))
    cat("    TLS score summary:\n")
    print(summary(tls_result$TLS.score))

    # Save
    saveRDS(se, file.path(OUT_DIR, paste0(samp$name_tag, "_seurat.rds")))
    saveRDS(tls_result, file.path(OUT_DIR, paste0(samp$name_tag, "_tls.rds")))

    # Plot
    pdf(file.path(OUT_DIR, paste0(samp$name_tag, "_TLS.pdf")), width = 20, height = 12)
    se$TLS_binary <- factor(ifelse(se$TLS_region == "TLS", "TLS", "Non-TLS"),
                            levels = c("TLS", "Non-TLS"))
    print(SpatialFeaturePlot(se, features = "TLS_score", pt.size.factor = 1.6) +
          ggplot2::labs(title = "TLS Score"))
    print(SpatialDimPlot(se, group.by = "TLS_binary", pt.size.factor = 1.6,
          cols = c("TLS" = "red", "Non-TLS" = "grey80")) +
          ggplot2::labs(title = paste0("TLS Regions (", n_tls, " spots)")))
    print(SpatialFeaturePlot(se, features = c("B_Plasma", "T_cell", "BT_codist", "LC50sig"),
          pt.size.factor = 1.6))
    print(SpatialDimPlot(se, group.by = "seurat_clusters", pt.size.factor = 1.6) +
          ggplot2::labs(title = "Seurat Clusters"))
    dev.off()

    results_list[[samp$name_tag]] <- data.frame(
      sample = samp$name_tag, total_spots = n_total,
      TLS_spots = n_tls, TLS_pct = round(100*n_tls/n_total, 2),
      max_TLS_score = round(max(tls_result$TLS.score, na.rm=TRUE), 4),
      max_B_Plasma = round(max(cell_scores$B_Plasma), 4),
      max_T_cell = round(max(cell_scores$T_cell), 4),
      max_BT_codist = round(max(b_t_codist), 4),
      stringsAsFactors = FALSE
    )
    cat("    Done!\n")

  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n")
  })

  unlink(tmp_dir, recursive = TRUE)
}

if (length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names = FALSE)
}
cat("\n====== Done ======\nResults:", OUT_DIR, "\n")
