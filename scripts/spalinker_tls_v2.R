# =============================================================================
# SpaLinker TLS Detection for GBM Visium Data (GSE237183)
# Full workflow: Cell annotation → TLS scoring → TLS spotting
# =============================================================================

library(Seurat)
library(SpaLinker)

# ---- Sample metadata (IDH-WT GBM only) ----
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

# Pilot: first 2
samples <- samples[1:2]

DATA_DIR <- "E:/GBM/GSE237183_RAW"
OUT_DIR  <- "E:/GBM/results/spalinker_tls_v2"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results_list <- list()

for (i in seq_along(samples)) {
  samp <- samples[[i]]
  cat(sprintf("\n========== [%d/%d] %s (%s) ==========\n",
              i, length(samples), samp$name_tag, samp$gsm))

  file_prefix <- file.path(DATA_DIR, paste0(samp$gsm, "_", samp$name_tag))
  h5_file <- paste0(file_prefix, "_filtered_feature_bc_matrix.h5")

  if (!file.exists(h5_file)) {
    cat("  SKIPPING: h5 file not found\n")
    next
  }

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
    # ---- Step 1: Load Visium data ----
    cat("  [1/6] Loading Visium data...\n")
    se <- Load10X_Spatial(data.dir = tmp_dir, assay = "Spatial",
                          filter.matrix = TRUE, slice = samp$name_tag)
    se$orig.ident <- samp$name_tag

    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)
    cat(sprintf("    Spots after QC: %d\n", ncol(se)))
    if (ncol(se) < 50) next

    # ---- Step 2: Preprocess (SCTransform + clustering) ----
    cat("  [2/6] Preprocessing (SCTransform, PCA, clustering)...\n")
    se <- SCTransform(se, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)
    se <- RunPCA(se, assay = "SCT", verbose = FALSE)
    se <- FindNeighbors(se, dims = 1:20, verbose = FALSE)
    se <- FindClusters(se, resolution = 0.6, verbose = FALSE)

    # ---- Step 3: Cell type annotation using PancancerCellAtlas ----
    cat("  [3/6] Computing cell type scores (PancancerCellAtlas 46 types)...\n")
    expr_counts <- GetAssayData(se, assay = "Spatial", slot = "counts")

    cell_scores <- GetGsetSigScore(
      expr     = expr_counts,
      category = "CuratedSig",
      types    = "Immune",
      subtype  = "PancancerCellAtlas",
      method   = "AddModuleScore",
      verbose  = FALSE
    )
    cell_mat <- cell_scores$Immune$PancancerCellAtlas
    cat(sprintf("    Cell type scores: %d spots x %d cell types\n", nrow(cell_mat), ncol(cell_mat)))

    # ---- Step 4: TLS gene signature scores ----
    cat("  [4/6] Computing TLS signature scores (8 TLS gene sets)...\n")
    tls_scores <- GetGsetSigScore(
      expr     = expr_counts,
      category = "CuratedSig",
      types    = "Immune",
      subtype  = "TLS",
      method   = "AddModuleScore",
      verbose  = FALSE
    )
    tls_mat <- tls_scores$Immune$TLS
    cat(sprintf("    TLS signature scores: %d spots x %d signatures\n", nrow(tls_mat), ncol(tls_mat)))

    # ---- Step 5: Create STFeature object ----
    cat("  [5/6] Adding spatial coordinates to meta.data...\n")
    # CreateStfObj requires x,y in se@meta.data - add from image coordinates
    coords <- se@images[[samp$name_tag]]@coordinates
    se$x <- coords$imagecol
    se$y <- coords$imagerow

    cat("      Creating STFeature object with cell annotations...\n")
    stf <- CreateStfObj(
      st        = se,
      cell.anno = cell_mat,
      init.fea  = c("Position", "CellAnno", "CellCodis"),
      assay     = "SCT",
      slot      = "data",
      norm      = TRUE,
      min.prop  = 0,
      verbose   = FALSE
    )

    # ---- Step 6: Run CalTLSfea ----
    cat("  [6/6] Running CalTLSfea for TLS detection...\n")

    # Use Seurat clusters for spatial neighborhood constraints
    clusters <- as.character(Idents(se))
    names(clusters) <- colnames(se)

    tls_result <- CalTLSfea(
      data       = tls_mat,           # TLS signature scores per spot
      st_pos     = stf@Position,      # Spatial coordinates
      cluster    = clusters,          # Seurat cluster for neighborhood
      cutoff     = 0.2,
      filt.dist  = 4,
      filt.spots = 2
    )

    # ---- Save results ----
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]

    n_tls <- sum(se$TLS_region == "TLS", na.rm = TRUE)
    n_total <- ncol(se)
    cat(sprintf("    TLS spots: %d / %d (%.1f%%)\n", n_tls, n_total, 100 * n_tls / n_total))

    # TLS score distribution
    score_summary <- summary(tls_result$TLS.score)
    cat("    TLS score summary:\n")
    print(score_summary)

    # Save Seurat object
    saveRDS(se, file.path(OUT_DIR, paste0(samp$name_tag, "_seurat.rds")))
    saveRDS(tls_result, file.path(OUT_DIR, paste0(samp$name_tag, "_tls_result.rds")))
    saveRDS(stf, file.path(OUT_DIR, paste0(samp$name_tag, "_stf.rds")))

    # Plot
    pdf(file.path(OUT_DIR, paste0(samp$name_tag, "_TLS.pdf")), width = 16, height = 5)

    if (n_tls > 0) {
      p1 <- SpatialFeaturePlot(se, features = "TLS_score", pt.size.factor = 1.6) +
        labs(title = paste0(samp$name_tag, " - TLS Score"))
      print(p1)
    }

    se$TLS_binary <- factor(ifelse(se$TLS_region == "TLS", "TLS", "Non-TLS"),
                            levels = c("TLS", "Non-TLS"))
    if (n_tls > 0) {
      p2 <- SpatialDimPlot(se, group.by = "TLS_binary", pt.size.factor = 1.6,
                           cols = c("TLS" = "red", "Non-TLS" = "grey80")) +
        labs(title = paste0(samp$name_tag, " - TLS Regions (", n_tls, " spots)"))
    } else {
      p2 <- SpatialDimPlot(se, group.by = "TLS_binary", pt.size.factor = 1.6,
                           cols = c("TLS" = "red", "Non-TLS" = "grey80")) +
        labs(title = paste0(samp$name_tag, " - No TLS detected"))
      # Also show top B cell and T cell scores
      if ("P_c0_CD8_KLRD1_GZMB_CCL4/5" %in% colnames(cell_mat)) {
        se$Bcell_Tcell <- cell_mat[colnames(se), "P_c0_CD8_KLRD1_GZMB_CCL4/5"]
        p3 <- SpatialFeaturePlot(se, features = "Bcell_Tcell", pt.size.factor = 1.6) +
          labs(title = paste0(samp$name_tag, " - B/plasma cell score"))
        print(p3)
      }
    }
    print(p2)

    # Cluster overview
    p4 <- SpatialDimPlot(se, group.by = "seurat_clusters", pt.size.factor = 1.6) +
      labs(title = paste0(samp$name_tag, " - Seurat Clusters"))
    print(p4)

    dev.off()

    # Summary
    results_list[[samp$name_tag]] <- data.frame(
      sample       = samp$name_tag,
      GSM          = samp$gsm,
      region       = samp$tumor_region,
      total_spots  = n_total,
      TLS_spots    = n_tls,
      TLS_pct      = round(100 * n_tls / n_total, 2),
      TLS_max_score = round(max(tls_result$TLS.score, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )

    cat("    Done!\n")

  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n")
  })

  unlink(tmp_dir, recursive = TRUE)
}

# ---- Final summary ----
if (length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names = FALSE)
}

cat("\n====== Done! ======\n")
cat("Results saved to:", OUT_DIR, "\n")
