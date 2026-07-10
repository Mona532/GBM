library(rhdf5)
library(Matrix)

root <- "E:/GBM"
h5_dir <- file.path(root, "spatial_data_visium", "spatial_data_visium", "anndata_consolidated")
tls_root <- file.path(root, "results", "tls_consolidated")
out_root <- file.path(root, "results")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

build_tls_components <- function(tls_df, radius = 2, min_spots = 5) {
  tls_only <- tls_df[tls_df[["TLS.region"]] == "TLS", , drop = FALSE]
  if (nrow(tls_only) < min_spots) return(NULL)

  xy <- cbind(
    as.numeric(tls_only[["array_col"]]),
    as.numeric(tls_only[["array_row"]]) * sqrt(3)
  )
  d <- as.matrix(dist(xy))
  adj <- d <= radius + 1e-8
  diag(adj) <- TRUE

  n <- nrow(adj)
  visited <- rep(FALSE, n)
  labels <- rep(NA_integer_, n)
  cid <- 0L

  for (i in seq_len(n)) {
    if (visited[i]) next
    cid <- cid + 1L
    queue <- i
    visited[i] <- TRUE
    members <- integer(0)

    while (length(queue) > 0) {
      v <- queue[1]
      queue <- queue[-1]
      members <- c(members, v)
      nei <- which(adj[v, ])
      new_nodes <- nei[!visited[nei]]
      if (length(new_nodes) > 0) {
        visited[new_nodes] <- TRUE
        queue <- c(queue, new_nodes)
      }
    }

    labels[members] <- cid - 1L
  }

  tls_only$component_id <- labels
  sizes <- table(tls_only$component_id)
  keep <- names(sizes[sizes >= min_spots])
  if (length(keep) == 0) return(NULL)
  tls_only[tls_only$component_id %in% keep, , drop = FALSE]
}

read_h5_counts <- function(path) {
  genes <- h5read(path, "var/_index")
  barcodes <- h5read(path, "obs/_index")
  xdata <- as.numeric(h5read(path, "X/data"))
  idx <- as.integer(h5read(path, "X/indices")) + 1L
  indptr <- as.integer(h5read(path, "X/indptr"))
  row_id <- rep(seq_len(length(barcodes)), diff(indptr))

  mat_obs_gene <- sparseMatrix(
    i = row_id,
    j = idx,
    x = xdata,
    dims = c(length(barcodes), length(genes))
  )
  rownames(mat_obs_gene) <- barcodes
  colnames(mat_obs_gene) <- genes

  q05 <- t(h5read(path, "obsm/c2l_ilc_q05"))
  ct <- h5read(path, "uns/c2l_ilc_cell_types")
  colnames(q05) <- ct
  rownames(q05) <- barcodes

  list(
    counts = t(mat_obs_gene),
    q05 = q05,
    genes = genes,
    barcodes = barcodes,
    celltypes = ct
  )
}

sample_ids <- intersect(
  sub("\\.h5ad$", "", basename(list.files(h5_dir, pattern = "\\.h5ad$", full.names = FALSE))),
  basename(list.dirs(tls_root, recursive = FALSE, full.names = FALSE))
)
sample_ids <- sort(sample_ids)

component_meta <- list()
component_spots <- list()
count_list <- list()
comp_list <- list()
all_ct <- NULL

for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[[idx]]
  tls_csv <- file.path(tls_root, sid, "tls_spot_scores_official_relaxed.csv")
  h5ad <- file.path(h5_dir, paste0(sid, ".h5ad"))
  if (!file.exists(tls_csv) || !file.exists(h5ad)) next

  tls_df <- read.csv(tls_csv, check.names = FALSE)
  comp_df <- build_tls_components(tls_df, radius = 2, min_spots = 5)
  if (is.null(comp_df)) next

  dat <- read_h5_counts(h5ad)
  if (is.null(all_ct)) all_ct <- dat$celltypes

  shared <- intersect(comp_df[["barcode"]], colnames(dat$counts))
  comp_df <- comp_df[comp_df[["barcode"]] %in% shared, , drop = FALSE]
  if (nrow(comp_df) == 0) next

  split_comp <- split(comp_df, comp_df[["component_id"]])
  kept <- 0L

  for (cid_name in names(split_comp)) {
    sub <- split_comp[[cid_name]]
    cid <- as.integer(cid_name)
    bc <- sub[["barcode"]]
    bc <- bc[bc %in% colnames(dat$counts)]
    if (length(bc) < 5) next

    unit_id <- paste0(sid, "__c", cid)
    count_list[[unit_id]] <- Matrix::rowSums(dat$counts[, bc, drop = FALSE])
    comp_list[[unit_id]] <- colMeans(dat$q05[bc, , drop = FALSE])
    kept <- kept + 1L

    component_meta[[unit_id]] <- data.frame(
      unit_id = unit_id,
      sample = sid,
      component_id = cid,
      n_spots = length(bc),
      tls_score_mean = mean(sub[["TLS.score"]], na.rm = TRUE),
      array_col_mean = mean(as.numeric(sub[["array_col"]]), na.rm = TRUE),
      array_row_mean = mean(as.numeric(sub[["array_row"]]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    component_spots[[unit_id]] <- data.frame(
      unit_id = unit_id,
      sample = sid,
      component_id = cid,
      barcode = bc,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("[%03d/%03d] %s: %d components kept\n", idx, length(sample_ids), sid, kept))
}

meta_df <- do.call(rbind, component_meta)
spots_df <- do.call(rbind, component_spots)
count_mat <- do.call(cbind, count_list)
comp_mat <- do.call(rbind, comp_list)
colnames(comp_mat) <- all_ct

saveRDS(count_mat, file.path(out_root, "tls_pseudobulk_counts_by_component.rds"))
saveRDS(comp_mat, file.path(out_root, "tls_pseudobulk_c2l_by_component.rds"))
write.csv(meta_df, file.path(out_root, "tls_pseudobulk_component_metadata.csv"), row.names = FALSE)
write.csv(spots_df, file.path(out_root, "tls_component_spot_map.csv"), row.names = FALSE)

summary_df <- data.frame(
  n_units = nrow(meta_df),
  n_samples = length(unique(meta_df$sample)),
  median_spots = median(meta_df$n_spots),
  mean_spots = mean(meta_df$n_spots),
  min_spots = min(meta_df$n_spots),
  max_spots = max(meta_df$n_spots)
)
write.csv(summary_df, file.path(out_root, "tls_pseudobulk_component_summary.csv"), row.names = FALSE)

cat("Saved TLS pseudobulk component outputs\n")
print(summary_df)
