# ============================================================
# joint_nmf.R — Joint NMF on composition + transcriptome
# ============================================================
# Concatenates C2L cell-type proportions with pseudobulk expression PCs,
# then runs NMF to discover programs defined by BOTH data modalities.
# Output: joint_nmf_rank5_weights.csv, joint_nmf_rank5_basis.csv
# ============================================================

suppressMessages({library(NMF); library(edgeR); library(matrixStats)})
root <- "E:/GBM/results"
comp   <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
wt     <- read.csv(file.path(root, "posdev_tau03_rank5_weights.csv"), check.names = FALSE)

common <- intersect(rownames(comp), colnames(counts))
comp   <- as.matrix(comp[common, , drop = FALSE]); storage.mode(comp) <- "numeric"
counts <- as.matrix(counts[, common, drop = FALSE])
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

# ---- Block 1: composition (proportions) ----
prop <- sweep(comp, 1, rowSums(comp) + 1e-8, "/")
prop <- prop[, colMeans(prop > 0) >= 0.05, drop = FALSE]   # drop rare cell types
prop <- scale(log1p(prop * 100))                            # log + scale per cell type
prop[!is.finite(prop)] <- 0

# ---- Block 2: transcriptome (top HVG → PCA) ----
logcpm <- cpm(calcNormFactors(DGEList(counts = counts)), log = TRUE, prior.count = 1)
rv     <- rowVars(logcpm)
top    <- names(sort(rv, decreasing = TRUE))[1:min(2000, length(rv))]
expr_z <- t(scale(t(logcpm[top, , drop = FALSE])))
expr_z[!is.finite(expr_z)] <- 0
pca    <- prcomp(t(expr_z), center = FALSE, scale. = FALSE)
expr_pc <- pca$x[, 1:min(20, ncol(pca$x)), drop = FALSE]

# ---- Concatenate blocks, scale each to equal total variance ----
X_comp <- scale(prop) / sqrt(ncol(prop))        # variance = 1 for block
X_expr <- scale(expr_pc) / sqrt(ncol(expr_pc))  # variance = 1 for block
X_joint <- cbind(X_comp, X_expr)                # components × features
X_joint[!is.finite(X_joint)] <- 0
X_nmf <- t(pmax(X_joint, 0))                    # features × components, non-negative

cat(sprintf("Joint NMF input: %d features × %d components (comp=%d + expr=%d)\n",
            nrow(X_nmf), ncol(X_nmf), ncol(prop), ncol(expr_pc)))

# ---- NMF rank 3-6 ----
set.seed(42)
for (k in 3:6) {
  fit <- nmf(X_nmf, rank = k, method = "brunet", nrun = 50, seed = 42)
  W <- basis(fit); H <- coef(fit)
  colnames(W) <- paste0("C", 1:k)
  rownames(H) <- paste0("C", 1:k); colnames(H) <- colnames(X_nmf)
  dom <- setNames(paste0("C", apply(H, 2, which.max)), colnames(X_nmf))

  cat(sprintf("\n--- rank %d ---\n", k))
  for (i in 1:k) {
    # composition signature (top cell types in W)
    w_comp <- W[colnames(prop), i]
    top_ct <- names(sort(w_comp, decreasing = TRUE))[1:5]

    # member components raw composition
    mem <- names(dom)[dom == paste0("C", i)]
    comp_raw <- colMeans(comp[mem, , drop = FALSE])
    top_raw <- names(sort(comp_raw, decreasing = TRUE))[1:5]

    cat(sprintf("C%d n=%d comp: %s\n       raw: %s\n",
                i, length(mem),
                paste(top_ct, collapse = ", "),
                paste(sprintf("%s(%.2f)", top_raw, comp_raw[top_raw]), collapse = ", ")))
  }

  wt <- as.data.frame(t(H))
  wt$unit_id <- rownames(wt)
  wt$dominant <- dom[rownames(wt)]
  wt <- wt[, c("unit_id", paste0("C", 1:k), "dominant")]
  write.csv(W,  file.path(root, paste0("joint_nmf_rank", k, "_basis.csv")))
  write.csv(wt, file.path(root, paste0("joint_nmf_rank", k, "_weights.csv")),
            row.names = FALSE)
}
cat("\nDone.\n")
