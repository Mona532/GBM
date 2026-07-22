library(NMF)

root <- "E:/GBM/results"
comp <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
meta <- read.csv(file.path(root, "tls_pseudobulk_component_metadata.csv"), check.names = FALSE)
meta$unit_id <- as.character(meta$unit_id)

comp <- as.matrix(comp)
storage.mode(comp) <- "numeric"
if (is.null(rownames(comp))) {
  stop("Component matrix is missing unit_id rownames; cannot align with metadata.")
}
shared <- intersect(meta$unit_id, rownames(comp))
if (length(shared) == 0) {
  stop("No overlapping unit_id between component matrix and metadata.")
}
meta <- meta[match(shared, meta$unit_id), , drop = FALSE]
comp <- comp[shared, , drop = FALSE]
stopifnot(identical(meta$unit_id, rownames(comp)))
mat <- t(comp)
mat <- pmax(mat, 0)

sample_n <- table(meta$sample)
sample_w <- 1 / as.numeric(sample_n[meta$sample])
names(sample_w) <- meta$unit_id
mat <- sweep(mat, 2, sample_w[colnames(mat)], "*")

normalize_factors <- function(W, H) {
  scale_vec <- colSums(W)
  scale_vec[scale_vec == 0] <- 1
  W_norm <- sweep(W, 2, scale_vec, "/")
  H_norm <- sweep(H, 1, scale_vec, "*")
  list(W = W_norm, H = H_norm, scale = scale_vec)
}

set.seed(42)
ranks <- 4:6
metrics <- list()
fit_list <- list()

for (k in ranks) {
  cat(sprintf("Running TLS component NMF rank %d...\n", k))
  fit <- nmf(mat, rank = k, method = "brunet", nrun = 50, seed = 42)
  fit_list[[as.character(k)]] <- fit
  s <- summary(fit)
  metrics[[as.character(k)]] <- data.frame(
    rank = k,
    cophenetic = s[["cophenetic"]],
    dispersion = s[["dispersion"]],
    sparseness_basis = s[["sparseness.basis"]],
    sparseness_coef = s[["sparseness.coef"]],
    silhouette_consensus = s[["silhouette.consensus"]],
    stringsAsFactors = FALSE
  )
}

metrics_df <- do.call(rbind, metrics)
write.csv(metrics_df, file.path(root, "tls_compnmf_rank_metrics.csv"), row.names = FALSE)

for (k in ranks) {
  fit <- fit_list[[as.character(k)]]
  W_raw <- basis(fit)
  H_raw <- coef(fit)
  norm <- normalize_factors(W_raw, H_raw)
  W <- norm$W
  H <- norm$H

  colnames(W) <- paste0("E", seq_len(ncol(W)))
  rownames(H) <- paste0("E", seq_len(nrow(H)))

  basis_df <- data.frame(
    cell_type = rownames(W),
    W,
    check.names = FALSE
  )
  write.csv(basis_df, file.path(root, paste0("tls_compnmf_rank", k, "_basis.csv")), row.names = FALSE)

  basis_raw_df <- data.frame(
    cell_type = rownames(W_raw),
    W_raw,
    check.names = FALSE
  )
  write.csv(basis_raw_df, file.path(root, paste0("tls_compnmf_rank", k, "_basis_raw.csv")), row.names = FALSE)

  weight_df <- as.data.frame(t(H))
  weight_df$unit_id <- rownames(weight_df)
  ecotype_cols <- paste0("E", seq_len(k))
  weight_df$dominant_ecotype <- ecotype_cols[max.col(weight_df[, ecotype_cols, drop = FALSE], ties.method = "first")]
  weight_df$dominant_weight <- apply(weight_df[, ecotype_cols, drop = FALSE], 1, max)
  weight_df$weight_total <- rowSums(weight_df[, ecotype_cols, drop = FALSE])
  weight_df <- merge(meta, weight_df, by = "unit_id", all.y = TRUE)
  write.csv(weight_df, file.path(root, paste0("tls_compnmf_rank", k, "_unit_weights.csv")), row.names = FALSE)

  raw_weight_df <- as.data.frame(t(H_raw))
  raw_weight_df$unit_id <- rownames(raw_weight_df)
  names(raw_weight_df)[seq_len(k)] <- ecotype_cols
  write.csv(raw_weight_df, file.path(root, paste0("tls_compnmf_rank", k, "_unit_weights_raw.csv")), row.names = FALSE)

  top_celltypes <- do.call(rbind, lapply(seq_len(k), function(i) {
    ord <- order(W[, i], decreasing = TRUE)
    topn <- head(ord, 8)
    data.frame(
      ecotype = paste0("E", i),
      cell_type = rownames(W)[topn],
      score = W[topn, i],
      rank = seq_along(topn),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(top_celltypes, file.path(root, paste0("tls_compnmf_rank", k, "_top_celltypes.csv")), row.names = FALSE)

  ecotype_meta <- do.call(rbind, lapply(seq_len(k), function(i) {
    sub <- weight_df[weight_df$dominant_ecotype == paste0("E", i), , drop = FALSE]
    if (nrow(sub) == 0) {
      return(data.frame(
        ecotype = paste0("E", i),
        n_units = 0,
        n_samples = 0,
        median_spots = NA,
        mean_tls_score = NA,
        dominant_weight_mean = NA,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      ecotype = paste0("E", i),
      n_units = nrow(sub),
      n_samples = length(unique(sub$sample)),
      median_spots = median(sub$n_spots),
      mean_tls_score = mean(sub$tls_score_mean, na.rm = TRUE),
      dominant_weight_mean = mean(sub$dominant_weight, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(ecotype_meta, file.path(root, paste0("tls_compnmf_rank", k, "_ecotype_summary.csv")), row.names = FALSE)
}

cat("Saved normalized TLS composition NMF outputs\n")
print(metrics_df)
