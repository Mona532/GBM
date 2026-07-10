# =============================================================================
# Reference-correlation deconvolution + SpaLinker TLS detection for GBM
# Uses scRNA-seq reference transcriptional profiles (not arbitrary genes)
# =============================================================================

library(Seurat)
library(SpaLinker)
library(ggplot2)
library(patchwork)

DATA_DIR <- "E:/GBM/GSE237183_RAW"
OUT_DIR  <- "E:/GBM/results/spalinker_tls_v4"
REF_PATH <- "E:/GBM/results/GBM_annotation/annotated_seurat.rds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Step 0: Build reference cell type profiles from scRNA-seq ----
cat("=== Building scRNA-seq reference profiles ===\n")
ref <- readRDS(REF_PATH)
ref <- NormalizeData(ref, verbose=FALSE)
ref_celltypes <- ref$llm_annotation_res0.5
cat("Cell types:", paste(names(table(ref_celltypes)), collapse=", "), "\n")

# Compute average expression per cell type (log-normalized)
Idents(ref) <- "llm_annotation_res0.5"
ref_avg <- AverageExpression(ref, assays="RNA", slot="data", verbose=FALSE)$RNA
cat(sprintf("Reference: %d genes x %d cell types\n", nrow(ref_avg), ncol(ref_avg)))

# For TLS: merge T cell subtypes into one "T cell" category
t_types <- c("T cell", "naive T cell", "CD8-positive, alpha-beta cytotoxic T cell")
ref_avg$T_cell <- rowMeans(ref_avg[, intersect(t_types, colnames(ref_avg)), drop=FALSE])
ref_avg <- ref_avg[, !colnames(ref_avg) %in% t_types]

# Ensure "B cell" column exists (SpaLinker naming requirement)
cat("Deconvolution cell types:", paste(colnames(ref_avg), collapse=", "), "\n\n")

# ---- Samples ----
samples <- list(
  list(gsm="GSM7596587", tag="mgh258",   region="unannotated"),
  list(gsm="GSM7596588", tag="zh881inf",  region="infiltrating"),
  list(gsm="GSM7596589", tag="zh881t1",   region="T1"),
  list(gsm="GSM7596590", tag="zh916bulk", region="bulk"),
  list(gsm="GSM7596591", tag="zh916inf",  region="infiltrating"),
  list(gsm="GSM7596592", tag="zh916t1",   region="T1"),
  list(gsm="GSM7596593", tag="zh1007inf", region="infiltrating"),
  list(gsm="GSM7596594", tag="zh1007nec", region="necrotic"),
  list(gsm="GSM7596595", tag="zh1019inf", region="infiltrating"),
  list(gsm="GSM7596596", tag="zh1019t1",  region="T1"),
  list(gsm="GSM7596597", tag="zh8811a",   region="bulk"),
  list(gsm="GSM7596598", tag="zh8811b",   region="bulk"),
  list(gsm="GSM7596599", tag="zh8812",    region="bulk")
)

# LC.50sig genes
lc50_genes <- CuratedSig.lt$Immune$TLS$LC.50sig

results_list <- list()

for(i in seq_along(samples)) {
  samp <- samples[[i]]
  cat(sprintf("=== [%d/%d] %s (%s) ===\n", i, length(samples), samp$tag, samp$region))

  # ---- Setup Visium directory ----
  file_prefix <- file.path(DATA_DIR, paste0(samp$gsm, "_", samp$tag))
  h5_file <- paste0(file_prefix, "_filtered_feature_bc_matrix.h5")
  if(!file.exists(h5_file)) { cat("  SKIP: no h5\n"); next }

  tmp_dir <- file.path(OUT_DIR, "tmp", samp$gsm)
  sp_dir  <- file.path(tmp_dir, "spatial")
  dir.create(sp_dir, recursive=TRUE, showWarnings=FALSE)
  file.copy(h5_file, file.path(tmp_dir, "filtered_feature_bc_matrix.h5"), overwrite=TRUE)
  for(ext in c("tissue_positions_list.csv.gz", "scalefactors_json.json.gz", "tissue_lowres_image.png.gz")) {
    src <- paste0(file_prefix, "_", ext)
    dsn <- file.path(sp_dir, sub("\\.gz$","",ext))
    if(file.exists(src)) R.utils::gunzip(src, destname=dsn, overwrite=TRUE, remove=FALSE)
  }

  tryCatch({
    # ---- Load spatial data ----
    se <- Load10X_Spatial(data.dir=tmp_dir, assay="Spatial", filter.matrix=TRUE, slice=samp$tag)
    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern="^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)
    cat(sprintf("  Spots: %d\n", ncol(se)))
    if(ncol(se) < 50) next

    # ---- Normalize spatial ----
    se <- SCTransform(se, assay="Spatial", verbose=FALSE, return.only.var.genes=FALSE)
    se <- RunPCA(se, assay="SCT", verbose=FALSE)
    se <- FindNeighbors(se, dims=1:20, verbose=FALSE)
    se <- FindClusters(se, resolution=0.6, verbose=FALSE)

    # ---- Reference-correlation deconvolution ----
    cat("  Deconvolution by reference correlation...\n")
    # Get SCT data for spatial
    sct_mat <- as.matrix(GetAssayData(se, assay="SCT", slot="data"))

    # Align genes
    common_genes <- intersect(rownames(sct_mat), rownames(ref_avg))
    cat(sprintf("  Common genes: %d\n", length(common_genes)))

    # Compute Spearman correlation per spot vs each reference cell type
    sct_sub <- sct_mat[common_genes, , drop=FALSE]
    ref_sub <- as.matrix(ref_avg[common_genes, , drop=FALSE])

    cell_props <- apply(sct_sub, 2, function(spot) {
      apply(ref_sub, 2, function(ct) cor(spot, ct, method="spearman"))
    })
    cell_props <- t(cell_props)  # spots x cell_types

    # Set negative correlations to 0, normalize each spot to sum=1
    cell_props[cell_props < 0] <- 0
    cell_props <- cell_props / rowSums(cell_props)
    cell_props[is.na(cell_props)] <- 0

    # Only keep B_cell and T_cell for SpaLinker cell.anno
    keep_types <- intersect(c("B cell", "T_cell"), colnames(cell_props))
    if(length(keep_types) < 2) {
      cat("  SKIP: missing B cell or T_cell in reference\n")
      next
    }
    cell_anno <- cell_props[, keep_types, drop=FALSE]
    cat(sprintf("  cell.anno: %d spots x %s\n", nrow(cell_anno), paste(colnames(cell_anno), collapse=", ")))

    # ---- Create STFeature with cell proportions ----
    cat("  Creating STFeature...\n")
    coords <- se@images[[samp$tag]]@coordinates
    se$x <- coords$imagecol
    se$y <- coords$imagerow

    stf <- CreateStfObj(st=se, cell.anno=cell_anno,
      init.fea=c("Position","CellAnno","CellCodis"),
      assay="SCT", slot="data", norm=TRUE, min.prop=0, verbose=FALSE)

    # ---- LC.50sig enrichment ----
    cat("  LC.50sig enrichment...\n")
    expr_counts <- GetAssayData(se, assay="Spatial", slot="counts")
    lc50_avail <- intersect(lc50_genes, rownames(expr_counts))
    tls_expr <- t(as.matrix(GetAssayData(se, assay="SCT", slot="data")[lc50_avail, , drop=FALSE]))
    lc50_score <- GsetScore(expr=expr_counts,
      geneset=list(LC50sig=lc50_avail), method="AddModuleScore",
      scale=FALSE, verbose=FALSE)$LC50sig

    # ---- Build TLS features ----
    # Get B cell_T cell co-distribution from STFeature
    bt_codist <- stf@CellCodis
    bt_col <- grep("B cell.*T_cell|T_cell.*B cell", colnames(bt_codist), value=TRUE)
    cat(sprintf("  B-T co-dist column: %s\n", paste(bt_col, collapse=", ")))
    if(length(bt_col) == 0) {
      cat("  SKIP: no B-T co-distribution column\n")
      next
    }

    tls_features <- data.frame(
      LC50sig    = lc50_score,
      BT_codist  = bt_codist[[bt_col[1]]],
      row.names  = colnames(se),
      check.names = FALSE
    )

    # ---- CalTLSfea ----
    cat("  Running CalTLSfea...\n")
    domains <- as.character(Idents(se))
    names(domains) <- colnames(se)
    st_pos <- stf@Position

    tls_result <- CalTLSfea(data=tls_features, st_pos=st_pos,
      cluster=domains, cutoff=0.2, filt.dist=4, filt.spots=2)

    # ---- Results ----
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]
    se$B_cell_prop <- cell_anno[colnames(se), "B cell"]
    se$T_cell_prop <- cell_anno[colnames(se), "T_cell"]
    se$BT_codist   <- bt_codist[[bt_col[1]]][colnames(se)]

    n_tls <- sum(se$TLS_region=="TLS", na.rm=TRUE)
    cat(sprintf("  TLS spots: %d / %d (%.1f%%)\n", n_tls, ncol(se), 100*n_tls/ncol(se)))
    cat("  TLS score summary:\n")
    print(summary(tls_result$TLS.score))

    # Save
    saveRDS(se,        file.path(OUT_DIR, paste0(samp$tag, "_seurat.rds")))
    saveRDS(tls_result, file.path(OUT_DIR, paste0(samp$tag, "_tls.rds")))

    # Plot (spot only)
    pdf(file.path(OUT_DIR, paste0(samp$tag, "_spot.pdf")), width=14, height=10)
    p1 <- SpatialFeaturePlot(se, features="TLS_score", pt.size.factor=2.0,
          image.alpha=0, stroke=0) +
          ggplot2::labs(title=paste(samp$tag, "- TLS Score"),
               subtitle=sprintf("max=%.3f, TLS spots=%d", max(se$TLS_score,na.rm=TRUE), n_tls)) +
          scale_fill_viridis_c(option="inferno")
    p2 <- SpatialFeaturePlot(se, features="B_cell_prop", pt.size.factor=2.0,
          image.alpha=0, stroke=0) +
          ggplot2::labs(title="B cell proportion") + scale_fill_viridis_c(option="magma")
    p3 <- SpatialFeaturePlot(se, features="T_cell_prop", pt.size.factor=2.0,
          image.alpha=0, stroke=0) +
          ggplot2::labs(title="T cell proportion") + scale_fill_viridis_c(option="magma")
    p4 <- SpatialFeaturePlot(se, features="BT_codist", pt.size.factor=2.0,
          image.alpha=0, stroke=0) +
          ggplot2::labs(title="B-T co-distribution") + scale_fill_viridis_c(option="plasma")
    print((p1 | p2) / (p3 | p4))
    dev.off()

    results_list[[samp$tag]] <- data.frame(
      sample=samp$tag, region=samp$region, spots=ncol(se),
      TLS_spots=n_tls, TLS_pct=round(100*n_tls/ncol(se),2),
      TLS_max=round(max(tls_result$TLS.score,na.rm=TRUE),4),
      B_max=round(max(se$B_cell_prop),4),
      T_max=round(max(se$T_cell_prop),4),
      BT_max=round(max(se$BT_codist),4),
      stringsAsFactors=FALSE
    )
    cat("  Done!\n\n")
  }, error=function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
  unlink(tmp_dir, recursive=TRUE)
}

# Summary
if(length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names=FALSE)
}
cat("\n===== Done =====\nResults:", OUT_DIR, "\n")
