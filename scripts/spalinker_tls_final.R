# =============================================================================
# SpaLinker TLS Detection — Complete Pipeline
# Uses pre-computed cell2location proportions from h5ad files
# Run: source("E:/GBM/scripts/spalinker_tls_final.R")
# =============================================================================

library(Seurat)
library(SpaLinker)
library(rhdf5)
library(Matrix)

# ---- Config ----
DATA_DIR <- "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"
OUT_DIR  <- "E:/GBM/results/spalinker_tls_final"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: Read cell2location cell state names from one file ----
cat("Reading cell state names...\n")
sample0 <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE)[1]
ft <- h5read(sample0, "var/feature_type")
cell_idx <- which(ft == "Cell state abundances")
vn <- h5read(sample0, "var/_index")
cell_states <- vn[cell_idx]
cat(sprintf("Found %d cell states:\n", length(cell_states)))
for(s in cell_states) cat(sprintf("  %s\n", s))

# ---- Step 2: Map cell states to SpaLinker naming ----
# SpaLinker needs: "B cell", "T cell" for CalCellCodis co-distribution
# Adjust these mappings based on actual cell state names
cell_type_map <- list()

# Auto-detect B cell and T cell states
b_states <- grep("B.cell|B_cell|Plasma|plasma|Bcell|B lymph", cell_states, value = TRUE, ignore.case = TRUE)
t_states <- grep("T.cell|T_cell|CD4|CD8|Tcell|T lymph", cell_states, value = TRUE, ignore.case = TRUE)

cat(sprintf("\nB/Plasma cell states: %d\n", length(b_states)))
for(s in b_states) cat(sprintf("  -> %s\n", s))
cat(sprintf("\nT cell states: %d\n", length(t_states)))
for(s in t_states) cat(sprintf("  -> %s\n", s))

# ---- Helper: read one h5ad sample ----
read_visium_h5ad <- function(fp) {
  # Read gene expression (subset)
  ft <- h5read(fp, "var/feature_type")
  gene_idx <- which(ft == "Gene Expression")
  cell_idx <- which(ft == "Cell state abundances")

  # Read gene expression matrix
  X_gene <- as(h5read(fp, "X", index = list(NULL, gene_idx)), "dgCMatrix")
  gene_names <- h5read(fp, "var/_index")[gene_idx]
  colnames(X_gene) <- gene_names

  # Read cell proportions
  cell_props <- as.matrix(h5read(fp, "X", index = list(NULL, cell_idx)))
  cell_names <- h5read(fp, "var/_index")[cell_idx]
  colnames(cell_props) <- cell_names

  # Read spatial coords
  spatial <- h5read(fp, "obsm/spatial")
  rownames(spatial) <- h5read(fp, "obs/_index")
  colnames(spatial) <- c("x", "y")

  # Read obs
  obs_idx <- h5read(fp, "obs/_index")
  sample_name <- h5read(fp, "obs/sample_name")
  if(is.factor(sample_name)) sample_name <- as.character(sample_name)

  list(
    counts = t(X_gene),
    cell_props = cell_props,
    spatial = spatial,
    barcodes = obs_idx,
    sample = sample_name[1]
  )
}

# ---- Step 3: Process each sample ----
h5ad_files <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE)
cat(sprintf("\nProcessing %d samples...\n", length(h5ad_files)))

# LC.50sig genes
lc50_genes <- CuratedSig.lt$Immune$TLS$LC.50sig

results_list <- list()

for(i in seq_along(h5ad_files)) {
  fp <- h5ad_files[[i]]
  cat(sprintf("\n[%d/%d] %s\n", i, length(h5ad_files), basename(fp)))

  tryCatch({
    # Read data
    dat <- read_visium_h5ad(fp)
    cat(sprintf("  Spots: %d, Genes: %d, Cell states: %d\n",
                nrow(dat$spatial), ncol(dat$counts), ncol(dat$cell_props)))

    # ---- Create Seurat object from gene expression ----
    se <- CreateSeuratObject(counts = dat$counts, assay = "Spatial",
                             meta.data = data.frame(row.names = dat$barcodes))
    se$x <- dat$spatial[, "x"]
    se$y <- dat$spatial[, "y"]

    # QC + normalize
    se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
    se <- subset(se, nFeature_Spatial > 200 & percent.mt < 30)
    cat(sprintf("  After QC: %d spots\n", ncol(se)))

    if(ncol(se) < 50) next

    se <- SCTransform(se, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)
    se <- RunPCA(se, assay = "SCT", verbose = FALSE)
    se <- FindNeighbors(se, dims = 1:20, verbose = FALSE)
    se <- FindClusters(se, resolution = 0.6, verbose = FALSE)

    # ---- Build cell.anno from cell2location proportions ----
    # Align cell_props with seurat cells
    cell_props <- dat$cell_props[colnames(se), , drop = FALSE]
    rownames(cell_props) <- colnames(se)

    # Keep only B and T cell columns for TLS
    keep_cols <- intersect(c(b_states, t_states), colnames(cell_props))
    if(length(keep_cols) < 2) {
      cat("  SKIP: missing B or T cell states\n")
      next
    }
    cell_anno <- cell_props[, keep_cols, drop = FALSE]
    cat(sprintf("  cell.anno: %s\n", paste(colnames(cell_anno), collapse = ", ")))

    # ---- Create STFeature ----
    stf <- CreateStfObj(st = se, cell.anno = cell_anno,
      init.fea = c("Position", "CellAnno", "CellCodis"),
      norm = TRUE, min.prop = 0, verbose = FALSE)

    # ---- LC.50sig enrichment ----
    expr <- GetAssayData(se, assay = "Spatial", slot = "counts")
    lc50_avail <- intersect(lc50_genes, rownames(expr))
    lc50_score <- GsetScore(expr = expr,
      geneset = list(LC50sig = lc50_avail),
      method = "AddModuleScore", scale = FALSE, verbose = FALSE)$LC50sig

    # ---- Get B-T co-distribution ----
    bt_codis <- stf@CellCodis
    bt_cols <- grep("B.*T|T.*B|Plasma.*T|T.*Plasma", colnames(bt_codis),
                    value = TRUE, ignore.case = TRUE)
    cat(sprintf("  B-T co-dist columns: %s\n", paste(bt_cols, collapse = ", ")))
    if(length(bt_cols) == 0) next
    bt_codist <- bt_codis[[bt_cols[1]]]
    names(bt_codist) <- rownames(bt_codis)

    # ---- Build TLS features ----
    tls_features <- data.frame(
      LC50sig    = lc50_score[colnames(se)],
      BT_codist  = bt_codist[colnames(se)],
      row.names  = colnames(se)
    )

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

    # Add results
    se$TLS_score  <- tls_result$TLS.score[colnames(se)]
    se$TLS_region <- tls_result$TLS.region[colnames(se)]

    n_tls <- sum(se$TLS_region == "TLS", na.rm = TRUE)
    cat(sprintf("  TLS spots: %d / %d (%.1f%%)\n", n_tls, ncol(se), 100*n_tls/ncol(se)))
    print(summary(tls_result$TLS.score))

    # Save
    tag <- gsub("\\.h5ad$", "", basename(fp))
    saveRDS(se, file.path(OUT_DIR, paste0(tag, "_seurat.rds")))
    saveRDS(tls_result, file.path(OUT_DIR, paste0(tag, "_tls.rds")))

    # Quick plot
    pdf(file.path(OUT_DIR, paste0(tag, "_TLS.pdf")), width = 14, height = 10)
    p1 <- SpatialFeaturePlot(se, features = "TLS_score", pt.size.factor = 2.0,
          image.alpha = 0, stroke = 0) +
          ggplot2::labs(title = paste(tag, "- TLS Score"),
            subtitle = sprintf("max=%.3f, TLS=%d", max(se$TLS_score, na.rm = TRUE), n_tls))
    print(p1)

    se$TLS_bin <- factor(ifelse(se$TLS_region == "TLS", "TLS", "Non-TLS"),
                         levels = c("TLS", "Non-TLS"))
    p2 <- SpatialDimPlot(se, group.by = "TLS_bin", pt.size.factor = 2.0,
          image.alpha = 0, stroke = 0,
          cols = c("TLS" = "red", "Non-TLS" = "grey80")) +
          ggplot2::labs(title = paste("TLS Regions (", n_tls, "spots)"))
    print(p2)
    dev.off()

    results_list[[tag]] <- data.frame(
      sample = tag, spots = ncol(se), TLS_spots = n_tls,
      TLS_pct = round(100*n_tls/ncol(se), 2),
      TLS_max = round(max(tls_result$TLS.score, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  })
  H5close()
}

# ---- Summary ----
if(length(results_list) > 0) {
  summ <- do.call(rbind, results_list)
  print(summ)
  write.csv(summ, file.path(OUT_DIR, "TLS_summary.csv"), row.names = FALSE)
}

cat("\n===== Done =====\nResults:", OUT_DIR, "\n")
