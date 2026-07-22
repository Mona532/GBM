# ============================================================
# assemble_seurat_rds.R — R-only: MEX → Seurat RDS + TLS metadata
# Run with R 4.3.2 (has Seurat v4 + SpaLinker)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root    <- "E:/GBM"
rds_dir <- file.path(root, "results", "visium_rds")
samples <- list.files(rds_dir, pattern = "^AT|^dryad|^GSE")

# TLS component → ecotype → spot
spot_map   <- read.csv(file.path(root, "results", "tls_component_spot_map.csv"),
                       stringsAsFactors = FALSE)
ecotype_wt <- read.csv(file.path(root, "results", "tls_compnmf_rank5_unit_weights.csv"),
                       stringsAsFactors = FALSE)
spot_map$ecotype <- ecotype_wt$dominant_ecotype[match(spot_map$unit_id, ecotype_wt$unit_id)]
spot_map$ecotype[is.na(spot_map$ecotype)] <- "none"

cat(sprintf("Found %d sample directories\n", length(samples)))
done <- 0; skip <- 0; err <- 0

for (sid in samples) {
  sample_dir <- file.path(rds_dir, sid)
  rds_file    <- file.path(sample_dir, paste0(sid, ".rds"))
  if (file.exists(rds_file)) { skip <- skip + 1; next }

  fm_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  sp_dir <- file.path(sample_dir, "spatial")
  if (!dir.exists(fm_dir) || !dir.exists(sp_dir)) next

  tryCatch({
    # Read MEX counts (matrix is spots×genes from Python, need genes×spots)
    mtx_path <- file.path(fm_dir, "matrix.mtx.gz")
    if (!file.exists(mtx_path)) next
    X <- readMM(mtx_path)
    # Transpose if needed: Read10X expects genes×barcodes
    if (ncol(X) > nrow(X) * 0.9) X <- t(X)   # more columns than rows → transpose
    bc <- readLines(file.path(fm_dir, "barcodes.tsv.gz"))
    ft <- readLines(file.path(fm_dir, "features.tsv.gz"))
    colnames(X) <- make.names(bc, unique = TRUE)
    gene_names <- sapply(strsplit(ft, "\t"), `[`, 1)
    rownames(X) <- make.names(gene_names, unique = TRUE)

    se <- CreateSeuratObject(counts = X, assay = "Spatial", project = sid)

    # Add hires image
    hires <- file.path(sp_dir, "tissue_hires_image.png")
    if (file.exists(hires)) {
      se@images[[sid]] <- Read10X_Image(
        image.dir  = sp_dir,
        image.name = "tissue_hires_image.png",
        filter.matrix = TRUE
      )
    }

    # Add TLS ecotype
    ss <- spot_map[spot_map$sample == sid, ]
    if (nrow(ss) > 0) {
      eco_vec <- setNames(ss$ecotype, ss$barcode)
      common  <- intersect(colnames(se), names(eco_vec))
      se$tls_ecotype <- "none"
      se$tls_ecotype[common] <- eco_vec[common]
    } else {
      se$tls_ecotype <- "none"
    }

    saveRDS(se, rds_file)
    done <- done + 1
    if (done %% 10 == 0) cat(sprintf("  %d/%d done\n", done, length(samples)))
  }, error = function(e) {
    err <<- err + 1
    cat(sprintf("%s: ERROR - %s\n", sid, e$message))
  })
}
cat(sprintf("Done: %d saved, %d skipped, %d errors\n", done, skip, err))
