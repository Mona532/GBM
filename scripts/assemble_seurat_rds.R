# ============================================================
# assemble_seurat_rds.R — Seurat RDS with TLS component metadata
# Reads spaceranger-layout dirs → Read10X + Read10X_Image → saveRDS
# ============================================================
suppressPackageStartupMessages({library(Seurat); library(Matrix)})

root    <- "E:/GBM"
rds_dir <- file.path(root, "results", "visium_rds")
samples <- list.files(rds_dir, pattern = "^AT|^dryad|^GSE")

# TLS component → ecotype → spot mapping
spot_map   <- read.csv(file.path(root, "results", "tls_component_spot_map.csv"),
                       stringsAsFactors = FALSE)
ecotype_wt <- read.csv(file.path(root, "results", "tls_compnmf_rank5_unit_weights.csv"),
                       stringsAsFactors = FALSE)
spot_map$ecotype <- ecotype_wt$dominant_ecotype[match(spot_map$unit_id, ecotype_wt$unit_id)]
spot_map$ecotype[is.na(spot_map$ecotype)] <- "none"

cat(sprintf("Found %d sample directories\n", length(samples)))

for (sid in samples) {
  sample_dir <- file.path(rds_dir, sid)
  rds_file    <- file.path(sample_dir, paste0(sid, ".rds"))
  if (file.exists(rds_file)) { cat(sprintf("%s: skip\n", sid)); next }

  fm_dir  <- file.path(sample_dir, "filtered_feature_bc_matrix")
  sp_dir  <- file.path(sample_dir, "spatial")
  if (!dir.exists(fm_dir) || !dir.exists(sp_dir)) next

  tryCatch({
    # Read counts (MEX format)
    counts <- Read10X(data.dir = fm_dir)
    se <- CreateSeuratObject(counts = counts, assay = "Spatial", project = sid)

    # Add spatial image (hires, skip if all-black)
    hires_path <- file.path(sp_dir, "tissue_hires_image.png")
    if (file.exists(hires_path)) {
      se@images[[sid]] <- Read10X_Image(
        image.dir  = sp_dir,
        image.name = "tissue_hires_image.png",
        filter.matrix = TRUE
      )
    }

    # Add TLS ecotype metadata
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
    cat(sprintf("%s: %d spots x %d genes, TLS=%d\n",
                sid, ncol(se), nrow(se), sum(se$tls_ecotype != "none")))
  }, error = function(e) {
    cat(sprintf("%s: ERROR - %s\n", sid, e$message))
  })
}
cat("Done\n")
