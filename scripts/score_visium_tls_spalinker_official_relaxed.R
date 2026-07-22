library(Matrix)
library(SpaLinker)
library(rhdf5)
library(Seurat)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript score_visium_tls_spalinker_official_relaxed.R <input_dir> <output_dir> [regex_pattern]")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
pattern <- if (length(args) >= 3) args[[3]] else "\\.h5ad$"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_h5ad_categorical <- function(file, group_path) {
  categories <- h5read(file, paste0(group_path, "/categories"))
  codes <- h5read(file, paste0(group_path, "/codes"))
  categories <- as.character(categories)
  codes <- as.integer(codes)
  out <- rep(NA_character_, length(codes))
  valid <- codes >= 0
  out[valid] <- categories[codes[valid] + 1L]
  out
}

read_h5ad_matrix <- function(file) {
  data <- h5read(file, "/X/data")
  indices <- h5read(file, "/X/indices")
  indptr <- h5read(file, "/X/indptr")
  n_obs <- length(h5read(file, "/obs/_index"))
  n_vars <- length(h5read(file, "/var/_index"))

  row_counts <- diff(indptr)
  i <- rep(seq_len(n_obs), row_counts)
  j <- as.integer(indices) + 1L
  sparseMatrix(
    i = i,
    j = j,
    x = as.numeric(data),
    dims = c(n_obs, n_vars)
  )
}

read_h5ad_minimal <- function(file) {
  obs_names <- as.character(h5read(file, "/obs/_index"))
  var_names <- as.character(h5read(file, "/var/_index"))
  feature_types <- read_h5ad_categorical(file, "/var/feature_types")
  coords <- t(h5read(file, "/obsm/spatial"))
  colnames(coords) <- c("x", "y")
  rownames(coords) <- obs_names

  obs <- data.frame(
    in_tissue = as.integer(h5read(file, "/obs/in_tissue")),
    array_row = as.integer(h5read(file, "/obs/array_row")),
    array_col = as.integer(h5read(file, "/obs/array_col")),
    sample = read_h5ad_categorical(file, "/obs/sample"),
    sample_name = read_h5ad_categorical(file, "/obs/sample_name"),
    row.names = obs_names,
    check.names = FALSE
  )

  x <- read_h5ad_matrix(file)
  rownames(x) <- obs_names
  colnames(x) <- var_names

  list(
    X = x,
    obs = obs,
    coords = coords,
    var_names = var_names,
    feature_types = feature_types
  )
}

sum_selected_features <- function(mat, names, pattern) {
  keep <- grepl(pattern, names, ignore.case = TRUE, perl = TRUE)
  if (!any(keep)) {
    return(list(values = rep(0, nrow(mat)), features = character(0)))
  }
  vals <- Matrix::rowSums(mat[, keep, drop = FALSE])
  list(values = as.numeric(vals), features = names[keep])
}

dominant_niche_cluster <- function(mat, feature_names) {
  if (ncol(mat) == 0) {
    return(rep("all", nrow(mat)))
  }
  idx <- max.col(as.matrix(mat), ties.method = "first")
  cluster_labels <- feature_names[idx]
  cluster_labels[cluster_labels == ""] <- "all"
  cluster_labels
}

get_tls_signature_scores <- function(expr) {
  data("CuratedSig.lt", package = "SpaLinker")
  geneset <- CuratedSig.lt[["Immune"]][["TLS"]]
  keep <- vapply(geneset, function(x) length(intersect(x, rownames(expr))) > 0, logical(1))
  geneset2 <- geneset[keep]
  temp <- CreateSeuratObject(counts = expr)
  n1 <- ncol(temp@meta.data)
  score_obj <- NULL
  avg_expr <- Matrix::rowMeans(GetAssayData(temp, assay = "RNA", slot = "counts"))
  candidate_bins <- unique(c(
    24L, 20L, 16L, 12L, 8L, 6L, 4L, 3L, 2L,
    min(24L, max(2L, length(unique(avg_expr)) - 1L))
  ))
  for (nbin in candidate_bins) {
    score_obj <- tryCatch(
      suppressWarnings(
        AddModuleScore(
          object = temp,
          assay = "RNA",
          slot = "counts",
          features = geneset2,
          name = names(geneset2),
          nbin = nbin,
          verbose = FALSE
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(score_obj)) {
      break
    }
  }
  if (is.null(score_obj)) {
    stop("Failed to calculate TLS signature scores with AddModuleScore")
  }
  n2 <- ncol(score_obj@meta.data)
  score <- score_obj@meta.data[, (n1 + 1):n2, drop = FALSE]
  colnames(score) <- names(geneset2)
  res <- matrix(NA_real_, ncol = length(geneset), nrow = nrow(score))
  rownames(res) <- rownames(score)
  colnames(res) <- names(geneset)
  res <- data.frame(res, check.names = FALSE)
  res[, colnames(score)] <- score
  res
}

score_one <- function(file, out_root) {
  sample_id <- sub("\\.h5ad$", "", basename(file))
  sample_dir <- file.path(out_root, sample_id)
  out_csv <- file.path(sample_dir, "tls_spot_scores_official_relaxed.csv")
  if (file.exists(out_csv)) {
    message("[skip] ", sample_id)
    spot_df <- read.csv(out_csv, check.names = FALSE)
    manifest_path <- file.path(sample_dir, "manifest.json")
    manifest <- if (file.exists(manifest_path)) jsonlite::fromJSON(manifest_path) else list()
    return(data.frame(
      sample_id = sample_id,
      n_spots = nrow(spot_df),
      n_tls_spots = sum(spot_df$TLS.region == "TLS"),
      tls_fraction = mean(spot_df$TLS.region == "TLS"),
      tls_score_mean = mean(spot_df$TLS.score, na.rm = TRUE),
      tls_score_median = median(spot_df$TLS.score, na.rm = TRUE),
      has_b_or_plasma = any(spot_df$Plasma_B_cells > 0),
      b_or_plasma_features = paste(manifest$b_or_plasma_features %||% character(0), collapse = "; "),
      t_features = paste(manifest$t_features %||% character(0), collapse = "; "),
      cluster_source = manifest$cluster_source %||% "dominant_spatial_niche",
      stringsAsFactors = FALSE
    ))
  }

  message("[score] ", sample_id)
  dat <- read_h5ad_minimal(file)

  gene_keep <- dat$feature_types == "Gene Expression"
  cell_keep <- dat$feature_types == "Cell state abundances"
  niche_keep <- dat$feature_types == "Spatial niche abundances"

  expr <- t(dat$X[, gene_keep, drop = FALSE])
  rownames(expr) <- make.unique(dat$var_names[gene_keep])
  colnames(expr) <- rownames(dat$obs)

  cell_mat <- dat$X[, cell_keep, drop = FALSE]
  colnames(cell_mat) <- dat$var_names[cell_keep]
  rownames(cell_mat) <- rownames(dat$obs)

  niche_mat <- dat$X[, niche_keep, drop = FALSE]
  colnames(niche_mat) <- dat$var_names[niche_keep]
  rownames(niche_mat) <- rownames(dat$obs)

  b_res <- sum_selected_features(cell_mat, colnames(cell_mat), "b cells|plasma")
  t_res <- sum_selected_features(cell_mat, colnames(cell_mat), "cd4|cd8|t cell")
  cluster <- dominant_niche_cluster(niche_mat, colnames(niche_mat))
  names(cluster) <- rownames(dat$obs)

  sig <- get_tls_signature_scores(expr)
  lc50 <- if ("LC.50sig" %in% colnames(sig)) {
    sig[, "LC.50sig"]
  } else {
    rep(0, ncol(expr))
  }
  names(lc50) <- colnames(expr)

  codis_in <- data.frame(
    "Plasma/B.cells" = b_res$values,
    "T.cells" = t_res$values,
    row.names = rownames(dat$obs),
    check.names = FALSE
  )
  codis <- CalCellCodis(codis_in, sort = TRUE)
  bt_col <- "Plasma/B.cells_T.cells"
  bt <- codis[, bt_col]

  tls_input <- data.frame(
    "Plasma/B.cells_T.cells" = bt,
    "LC.50sig" = lc50[rownames(dat$obs)],
    row.names = rownames(dat$obs),
    check.names = FALSE
  )
  st_pos <- data.frame(dat$coords[, c("x", "y"), drop = FALSE], check.names = FALSE)

  tls <- CalTLSfea(
    data = tls_input,
    st_pos = st_pos,
    cluster = cluster,
    cutoff = 0.2,
    filt.dist = 4,
    filt.spots = 2,
    r.dist = 2,
    method = "weighted",
    layer.method = "mean",
    adjust = 2,
    verbose = FALSE
  )

  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  spot_df <- data.frame(
    barcode = rownames(dat$obs),
    sample_id = sample_id,
    x = dat$coords[, "x"],
    y = dat$coords[, "y"],
    in_tissue = dat$obs$in_tissue,
    array_row = dat$obs$array_row,
    array_col = dat$obs$array_col,
    dominant_niche_cluster = cluster[rownames(dat$obs)],
    Plasma_B_cells = codis_in[, "Plasma/B.cells"],
    T_cells = codis_in[, "T.cells"],
    Plasma_B.cells_T.cells = bt[rownames(dat$obs)],
    LC.50sig = lc50[rownames(dat$obs)],
    TLS.score = tls$TLS.score[rownames(dat$obs)],
    TLS.region = tls$TLS.region[rownames(dat$obs)],
    check.names = FALSE
  )
  write.csv(spot_df, file.path(sample_dir, "tls_spot_scores_official_relaxed.csv"), row.names = FALSE)

  summary_row <- data.frame(
    sample_id = sample_id,
    n_spots = nrow(spot_df),
    n_tls_spots = sum(spot_df$TLS.region == "TLS"),
    tls_fraction = mean(spot_df$TLS.region == "TLS"),
    tls_score_mean = mean(spot_df$TLS.score, na.rm = TRUE),
    tls_score_median = median(spot_df$TLS.score, na.rm = TRUE),
    has_b_or_plasma = any(codis_in[, "Plasma/B.cells"] > 0),
    b_or_plasma_features = paste(b_res$features, collapse = "; "),
    t_features = paste(t_res$features, collapse = "; "),
    cluster_source = "dominant_spatial_niche",
    stringsAsFactors = FALSE
  )

  meta <- list(
    sample_id = sample_id,
    source_h5ad = file,
    source_logic = "SpaLinker::CalTLSfea official implementation with relaxed Plasma/B.cells input",
    official_inputs = c("Plasma/B.cells_T.cells", "LC.50sig"),
    b_or_plasma_features = b_res$features,
    t_features = t_res$features,
    cluster_source = "dominant Spatial niche abundances",
    cutoff = 0.2,
    filt.dist = 4,
    filt.spots = 2,
    regfea_args = list(r.dist = 2, method = "weighted", layer.method = "mean", adjust = 2)
  )
  write(jsonlite::toJSON(meta, pretty = TRUE, auto_unbox = TRUE), file.path(sample_dir, "manifest.json"))

  summary_row
}

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
files <- sort(files)
if (!length(files)) {
  stop("No .h5ad files found under input_dir")
}

summary_df <- do.call(rbind, lapply(files, score_one, out_root = output_dir))
write.csv(summary_df, file.path(output_dir, "tls_official_relaxed_summary.csv"), row.names = FALSE)
message("[done] ", nrow(summary_df), " samples written to ", output_dir)
