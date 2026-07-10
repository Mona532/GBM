# ============================================================
# 01_posdev_nmf.R — TLS component compositional program discovery
# ============================================================
# Input:  tls_pseudobulk_c2l_by_component.rds (C2L abundance, component x cell_type)
#         tls_pseudobulk_counts_by_component.rds (pseudobulk counts, gene x component)
# Output: posdev_tau03_rank5_weights.csv (per-component NMF factor weights + dominant)
#         posdev_tau03_rank5_basis.csv  (cell_type x factor basis matrix)
# Method: positive-deviation NMF (z-score → pmax(z-τ,0) → NMF)
#         τ = 0.3 filters weak enrichments, retains only robust positive deviations
# ============================================================

suppressMessages(library(NMF))

root  <- "E:/GBM/results"
comp  <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
# comp is component × cell_type; convert to numeric, drop empty rows/cols
comp  <- as.matrix(comp); storage.mode(comp) <- "numeric"
comp  <- comp[rowSums(comp) > 0, colSums(comp) > 0, drop = FALSE]

# ---- positive-deviation transformation ----
# 1. per-component proportion (removes density gradient, non-negative)
prop <- sweep(comp, 1, rowSums(comp) + 1e-8, "/")
# 2. drop cell types present in <5% of components (too rare to be informative)
prop <- prop[, colMeans(prop > 0) >= 0.05, drop = FALSE]
# 3. log-transform to compress high-abundance cell types
x <- log1p(prop * 100)
# 4. per-cell-type z-score across all components (standardizes scale)
z <- scale(x); z[!is.finite(z)] <- 0
# 5. threshold: only keep values > τ standard deviations above mean
tau <- 0.3
x_pos <- pmax(z - tau, 0)
# 6. drop cell types with zero positive entries after thresholding
x_pos <- x_pos[, colSums(x_pos) > 0, drop = FALSE]
# 7. transpose to cell_type × component for NMF
x_nmf <- t(x_pos)

cat(sprintf("NMF input: %d cell_types × %d components (tau=%.1f, nonzero=%.0f%%)\n",
            nrow(x_nmf), ncol(x_nmf), tau, 100 * mean(x_nmf > 0)))

# ---- NMF decomposition, rank 4-6 ----
set.seed(42)
for (k in c(4, 5, 6)) {
  fit <- nmf(x_nmf, rank = k, method = "brunet", nrun = 50, seed = 42)
  W <- basis(fit); H <- coef(fit)
  colnames(W) <- paste0("C", 1:k)
  rownames(H) <- paste0("C", 1:k)
  colnames(H) <- colnames(x_nmf)

  # dominant factor per component
  dom <- apply(H, 2, which.max)
  dominant <- setNames(paste0("C", dom), colnames(x_nmf))

  # full H weights (not just dominant) for downstream mixture analysis
  wt <- as.data.frame(t(H))
  wt$unit_id <- rownames(wt)
  wt$dominant <- dominant[rownames(wt)]
  wt <- wt[, c("unit_id", paste0("C", 1:k), "dominant")]

  write.csv(W,  file.path(root, paste0("posdev_tau03_rank", k, "_basis.csv")))
  write.csv(wt, file.path(root, paste0("posdev_tau03_rank", k, "_weights.csv")),
            row.names = FALSE)

  # report basis (top cell types per factor) and raw composition (for naming)
  cat(sprintf("\n--- rank %d ---\n", k))
  for (i in 1:k) {
    top_basis <- names(sort(W[, i], decreasing = TRUE))[1:5]
    mem  <- names(dom)[dom == i]
    comp_raw <- colMeans(comp[mem, , drop = FALSE])
    top_comp  <- names(sort(comp_raw, decreasing = TRUE))[1:5]
    cat(sprintf("C%d n=%d basis: %s\n       comp: %s\n",
                i, length(mem),
                paste(top_basis, collapse = ", "),
                paste(sprintf("%s(%.2f)", top_comp, comp_raw[top_comp]),
                      collapse = ", ")))
  }
}
cat("\nDone.\n")
