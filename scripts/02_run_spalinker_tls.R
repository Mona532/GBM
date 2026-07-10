# =============================================================================
# SpaLinker TLS Detection — processes all extracted h5ad samples
# Run: Rscript E:/GBM/scripts/02_run_spalinker_tls.R
# =============================================================================

library(Seurat)
library(SpaLinker)
library(Matrix)

INPUT_DIR <- "E:/GBM/results/spalinker_input"
OUT_DIR   <- "E:/GBM/results/spalinker_tls_final"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Read cell state mapping ----
info_file <- file.path(INPUT_DIR, "cell_state_info.txt")
lines <- readLines(info_file)
b_start <- grep("B_PLASMA_STATES:", lines) + 1
t_start <- grep("T_CELL_STATES:", lines) + 1

b_states <- trimws(lines[b_start:(t_start-2)])
t_states <- trimws(lines[t_start:length(lines)])
b_states <- b_states[b_states != ""]
t_states <- t_states[t_states != ""]

cat("B/Plasma cell states:\n")
for(s in b_states) cat(sprintf("  %s\n", s))
cat("T cell states:\n")
for(s in t_states) cat(sprintf("  %s\n", s))

# ---- LC.50sig genes ----
lc50_genes <- CuratedSig.lt$Immune$TLS$LC.50sig

# ---- Process all samples ----
sample_dirs <- list.dirs(INPUT_DIR, recursive = FALSE, full.names = TRUE)
cat(sprintf("\nProcessing %d samples...\n", length(sample_dirs)))

results_list <- list()

for(i in seq_along(sample_dirs)) {
  sdir <- sample_dirs[[i]]
  tag <- basename(sdir)
  cat(sprintf("\n[%d/%d] %s\n", i, length(sample_dirs), tag))

  # Check for required files
  mtx_file  <- file.path(sdir, "counts.mtx")
  gene_file <- file.path(sdir, "genes.tsv")
  bc_file   <- file.path(sdir, "barcodes.tsv")
  sp_file   <- file.path(sdir, "spatial.tsv")
  cp_file   <- file.path(sdir, "cell_props.tsv")

  if(!all(file.exists(mtx_file, gene_file, bc_file, sp_file, cp_file))) {
    cat("  SKIP: missing files\n")
    next
  }

  tryCatch({
    # ---- Load data ----
    counts <- readMM(mtx_file)
    genes  <- readLines(gene_file)
    barcodes <- readLines(bc_file)
    spatial  <- read.table(sp_file, header = TRUE)
    cell_props_raw <- as.matrix(read.table(cp_file))

    # Get cell state names from first sample's info
    all_cell_states <- trimws(lines[2:(b_start-2)])
    all_cell_states <- all_cell_states[all_cell_states != ""]
    colnames(cell_props_raw) <- all_cell_states

    rownames(counts) <- genes
    colnames(counts)  <- barcodes
    rownames(spatial) <- barcodes
    rownames(cell_props_raw) <- barcodes

    cat(sprintf("  Spots: %d, Genes: %d, Cell states: %d\n",
                length(barcodes), length(genes), ncol(cell_props_raw)))

    # ---- Build Seurat ----
    se <- CreateSeuratObject(counts = counts, assay = "Spatial",
                             meta.data = data.frame(row.names = barcodes))
    se$x <- spatial[barcodes, "x"]
    se$y <- spatial[barcodes, "y"]

    # QC
    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)
    cat(sprintf("  After QC: %d spots\n", ncol(se)))
    if(ncol(se) < 50) next

    # Normalize
    se <- SCTransform(se, assay = "Spatial", verbose = FALSE,
                      return.only.var.genes = FALSE)
    se <- RunPCA(se, assay = "SCT", verbose = FALSE)
    se <- FindNeighbors(se, dims = 1:20, verbose = FALSE)
    se <- FindClusters(se, resolution = 0.6, verbose = FALSE)

    # ---- Build cell.anno: B + T cell proportions ----
    keep_types <- intersect(c(b_states, t_states), colnames(cell_props_raw))
    if(length(keep_types) < 2) { cat("  SKIP: insufficient B/T types\n"); next }

    cell_anno <- cell_props_raw[colnames(se), keep_types, drop = FALSE]
    cell_anno <- cell_anno / rowSums(cell_anno)
    cell_anno[is.na(cell_anno)] <- 0

    cat(sprintf("  cell.anno types: %s\n", paste(colnames(cell_anno), collapse = ", ")))

    # ---- Create STFeature ----
    stf <- CreateStfObj(st = se, cell.anno = cell_anno,
      init.fea = c("Position", "CellAnno", "CellCodis"),
      norm = TRUE, min.prop = 0, verbose = FALSE)

    # Get B-T co-distribution columns
    bt_cols <- colnames(stf@CellCodis)
    bt_match <- apply(combn(colnames(cell_anno), 2), 2, function(pair) {
      paste(sort(pair), collapse = "_")
    })
    bt_codis <- stf@CellCodis
    cat(sprintf("  Co-dist columns: %d\n", ncol(bt_codis)))

    # ---- LC.50sig ----
    expr <- GetAssayData(se, assay = "Spatial", slot = "counts")
    lc50_avail <- intersect(lc50_genes, rownames(expr))
    lc50_score <- GsetScore(expr = expr,
      geneset = list(LC50sig = lc50_avail),
      method = "AddModuleScore", scale = FALSE, verbose = FALSE)$LC50sig

    # ---- Build TLS features ----
    # Use all pairwise co-distribution cols with LC.50sig
    tls_features <- bt_codis[colnames(se), , drop = FALSE]
    tls_features$LC50sig <- lc50_score[colnames(se)]

    cat(sprintf("  TLS features: %d columns\n", ncol(tls_features)))

    # ---- CalTLSfea ----
    domains <- as.character(Idents(se))
    names(domains) <- colnames(se)

    tls_result <- CalTLSfea(
      data       = tls_features,
      st_pos     = stf@Position,
      cluster    = domains,
      cutoff     = 0.2,
      filt.dist  = 4,
      filt.spots = 2
    )

    # Results
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]
    n_tls <- sum(se$TLS_region == "TLS", na.rm = TRUE)
    cat(sprintf("  TLS spots: %d / %d (%.1f%%)\n", n_tls, ncol(se), 100*n_tls/ncol(se)))
    print(summary(tls_result$TLS.score))

    # Save
    out_rds <- file.path(OUT_DIR, paste0(tag, "_seurat.rds"))
    saveRDS(se, out_rds)
    saveRDS(tls_result, file.path(OUT_DIR, paste0(tag, "_tls.rds")))

    # Quick spot-only plot
    pdf(file.path(OUT_DIR, paste0(tag, "_TLS.pdf")), width = 14, height = 5)
    p1 <- SpatialFeaturePlot(se, features = "TLS_score", pt.size.factor = 2.0,
          image.alpha = 0, stroke = 0) +
          ggplot2::labs(title = paste(tag, "- TLS Score"),
            subtitle = sprintf("TLS spots=%d, max=%.3f", n_tls, max(se$TLS_score,na.rm=TRUE)))
    print(p1)

    se$TLS_bin <- factor(ifelse(se$TLS_region=="TLS","TLS","Non-TLS"),
                         levels=c("TLS","Non-TLS"))
    cols <- if(n_tls > 0) c("TLS"="red","Non-TLS"="grey80") else c("TLS"="red","Non-TLS"="grey80")
    p2 <- SpatialDimPlot(se, group.by="TLS_bin", pt.size.factor=2.0,
          image.alpha=0, stroke=0, cols=cols) +
          ggplot2::labs(title = paste("TLS Regions:", n_tls, "spots"))
    print(p2)
    dev.off()

    results_list[[tag]] <- data.frame(
      sample=tag, spots=ncol(se), TLS_spots=n_tls,
      TLS_pct=round(100*n_tls/ncol(se),2),
      TLS_max=round(max(tls_result$TLS.score,na.rm=TRUE),4),
      stringsAsFactors=FALSE
    )

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  })
}

# ---- Summary ----
if(length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names = FALSE)
  cat(sprintf("\nProcessed %d/%d samples. %d with TLS>0.\n",
              nrow(summ), length(sample_dirs), sum(summ$TLS_spots > 0)))
}

cat("\n===== Done =====\n")
cat("Results:", OUT_DIR, "\n")
