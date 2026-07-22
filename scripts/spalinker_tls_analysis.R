# =============================================================================
# SpaLinker TLS Spot Identification for GBM Visium Data (GSE237183)
# =============================================================================

library(Seurat)
library(SpaLinker)

# ---- Sample metadata (IDH-WT GBM only) ----
# name_tag must match the file prefix in the RAW directory
samples <- list(
  list(gsm = "GSM7596587", name_tag = "mgh258",   tumor_region = "unannotated",  mgmt = "methylated"),
  list(gsm = "GSM7596588", name_tag = "zh881inf",  tumor_region = "infiltrating", mgmt = "partial"),
  list(gsm = "GSM7596589", name_tag = "zh881t1",   tumor_region = "T1",           mgmt = "partial"),
  list(gsm = "GSM7596590", name_tag = "zh916bulk", tumor_region = "bulk",         mgmt = "methylated"),
  list(gsm = "GSM7596591", name_tag = "zh916inf",  tumor_region = "infiltrating", mgmt = "methylated"),
  list(gsm = "GSM7596592", name_tag = "zh916t1",   tumor_region = "T1",           mgmt = "methylated"),
  list(gsm = "GSM7596593", name_tag = "zh1007inf", tumor_region = "infiltrating", mgmt = "NA"),
  list(gsm = "GSM7596594", name_tag = "zh1007nec", tumor_region = "necrotic",     mgmt = "NA"),
  list(gsm = "GSM7596595", name_tag = "zh1019inf", tumor_region = "infiltrating", mgmt = "NA"),
  list(gsm = "GSM7596596", name_tag = "zh1019t1",  tumor_region = "T1",           mgmt = "NA"),
  list(gsm = "GSM7596597", name_tag = "zh8811a",   tumor_region = "bulk",         mgmt = "partial"),
  list(gsm = "GSM7596598", name_tag = "zh8811b",   tumor_region = "bulk",         mgmt = "partial"),
  list(gsm = "GSM7596599", name_tag = "zh8812",    tumor_region = "bulk",         mgmt = "partial")
)

# Pilot: only first 2
samples <- samples[1:2]

DATA_DIR <- "E:/GBM/GSE237183_RAW"
OUT_DIR  <- "E:/GBM/results/spalinker_tls"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- TLS signature genes ----
tls_sigs <- CuratedSig.lt$Immune$TLS
cat("Available TLS signatures:", names(tls_sigs), "\n")
cat("Total unique TLS genes:", length(unique(unlist(tls_sigs))), "\n")

# ---- Process each sample ----
results_list <- list()

for (i in seq_along(samples)) {
  samp <- samples[[i]]
  cat(sprintf("\n========== [%d/%d] %s (%s) ==========\n",
              i, length(samples), samp$name_tag, samp$gsm))

  file_prefix <- file.path(DATA_DIR, paste0(samp$gsm, "_", samp$name_tag))
  h5_file <- paste0(file_prefix, "_filtered_feature_bc_matrix.h5")

  if (!file.exists(h5_file)) {
    cat("  SKIPPING: h5 file not found at", h5_file, "\n")
    next
  }

  # Build proper Visium directory structure for Load10X_Spatial
  tmp_dir <- file.path(OUT_DIR, "tmp", samp$gsm)
  spatial_subdir <- file.path(tmp_dir, "spatial")
  dir.create(spatial_subdir, recursive = TRUE, showWarnings = FALSE)

  file.copy(h5_file, file.path(tmp_dir, "filtered_feature_bc_matrix.h5"),
            overwrite = TRUE)

  tissue_pos_gz <- paste0(file_prefix, "_tissue_positions_list.csv.gz")
  if (file.exists(tissue_pos_gz)) {
    R.utils::gunzip(tissue_pos_gz,
      destname = file.path(spatial_subdir, "tissue_positions_list.csv"),
      overwrite = TRUE, remove = FALSE)
  }

  scalefactors_gz <- paste0(file_prefix, "_scalefactors_json.json.gz")
  if (file.exists(scalefactors_gz)) {
    R.utils::gunzip(scalefactors_gz,
      destname = file.path(spatial_subdir, "scalefactors_json.json"),
      overwrite = TRUE, remove = FALSE)
  }

  tissue_img_gz <- paste0(file_prefix, "_tissue_lowres_image.png.gz")
  if (file.exists(tissue_img_gz)) {
    R.utils::gunzip(tissue_img_gz,
      destname = file.path(spatial_subdir, "tissue_lowres_image.png"),
      overwrite = TRUE, remove = FALSE)
  }

  tryCatch({
    cat("  Loading Visium data...\n")
    se <- Load10X_Spatial(data.dir = tmp_dir, assay = "Spatial",
                          filter.matrix = TRUE, slice = samp$name_tag)
    se$orig.ident <- samp$name_tag
    se$tumor_region <- samp$tumor_region
    se$mgmt_status <- samp$mgmt

    # Basic QC
    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)

    cat(sprintf("  Spots after QC: %d\n", ncol(se)))

    if (ncol(se) < 50) {
      cat("  SKIPPING: too few spots after QC\n")
      next
    }

    # Preprocess with SCTransform
    cat("  Running SCTransform...\n")
    se <- SCTransform(se, assay = "Spatial", verbose = FALSE)

    # Ensure SCT assay exists
    if (!"SCT" %in% names(se@assays)) {
      cat("  ERROR: SCT assay not created\n")
      next
    }

    # Clustering
    se <- RunPCA(se, assay = "SCT", verbose = FALSE)
    se <- FindNeighbors(se, dims = 1:20, verbose = FALSE)
    se <- FindClusters(se, resolution = 0.6, verbose = FALSE)

    # Get SCT data using GetAssayData (safe for S4 objects)
    sct_data <- GetAssayData(se, assay = "SCT", slot = "data")
    cat(sprintf("  SCT data dimensions: %d genes x %d spots\n",
                nrow(sct_data), ncol(sct_data)))

    # Collect TLS signature genes
    tls_all_genes <- unique(unlist(tls_sigs))
    tls_genes_available <- intersect(tls_all_genes, rownames(sct_data))
    cat(sprintf("  TLS signature genes available: %d / %d\n",
                length(tls_genes_available), length(tls_all_genes)))

    if (length(tls_genes_available) < 5) {
      cat("  SKIPPING: too few TLS genes available\n")
      next
    }

    # Expression matrix: spots x genes (required by CalTLSfea / RegFeaEnrich)
    expr_mat <- t(as.matrix(sct_data[tls_genes_available, , drop = FALSE]))

    # Spatial positions (image pixel coordinates)
    st_pos <- se@images[[samp$name_tag]]@coordinates[, c("imagecol", "imagerow")]
    colnames(st_pos) <- c("x", "y")
    rownames(st_pos) <- colnames(se)

    # Cluster assignments
    clusters <- as.character(Idents(se))
    names(clusters) <- colnames(se)

    # ----- Run CalTLSfea -----
    cat("  Running CalTLSfea...\n")
    tls_result <- CalTLSfea(
      data       = expr_mat,
      st_pos     = st_pos,
      cluster    = clusters,
      cutoff     = 0.2,
      filt.dist  = 4,
      filt.spots = 2
    )

    # Add results to Seurat object
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]

    n_tls <- sum(se$TLS_region == "TLS", na.rm = TRUE)
    n_total <- ncol(se)
    cat(sprintf("  TLS spots: %d / %d (%.1f%%)\n",
                n_tls, n_total, 100 * n_tls / n_total))
    cat(sprintf("  TLS score summary: min=%.4f, q25=%.4f, median=%.4f, q75=%.4f, max=%.4f\n",
                quantile(tls_result$TLS.score, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)))

    # Save spatial PDF using labs() instead of ggtitle
    pdf(file.path(OUT_DIR, paste0(samp$name_tag, "_TLS.pdf")), width = 16, height = 5)

    p1 <- SpatialFeaturePlot(se, features = "TLS_score", pt.size.factor = 1.6) +
      labs(title = paste0(samp$name_tag, " - TLS Score"))

    se$TLS_binary <- ifelse(se$TLS_region == "TLS", "TLS", "Non-TLS")
    p2 <- SpatialDimPlot(se, group.by = "TLS_binary", pt.size.factor = 1.6,
                         cols = c("TLS" = "red", "Non-TLS" = "grey80")) +
      labs(title = paste0(samp$name_tag, " - TLS Regions (", n_tls, " spots)"))

    p3 <- SpatialDimPlot(se, group.by = "seurat_clusters", pt.size.factor = 1.6) +
      labs(title = paste0(samp$name_tag, " - Clusters"))

    print(p1)
    print(p2)
    print(p3)
    dev.off()

    # Save results as RDS
    saveRDS(list(sample = samp, se = se, tls_result = tls_result,
                 tls_genes = tls_genes_available),
            file.path(OUT_DIR, paste0(samp$name_tag, "_TLS_result.rds")))

    results_list[[samp$name_tag]] <- data.frame(
      sample      = samp$name_tag,
      GSM         = samp$gsm,
      region      = samp$tumor_region,
      mgmt        = samp$mgmt,
      total_spots = n_total,
      TLS_spots   = n_tls,
      TLS_pct     = round(100 * n_tls / n_total, 2),
      n_TLS_genes = length(tls_genes_available),
      stringsAsFactors = FALSE
    )

    cat("  Done!\n")

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })

  # Clean up temp files
  unlink(tmp_dir, recursive = TRUE)
}

# ---- Summary ----
if (length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names = FALSE)
}

cat("\n====== Done! ======\n")
cat("Results saved to:", OUT_DIR, "\n")
