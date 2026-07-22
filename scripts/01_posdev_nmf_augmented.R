suppressMessages(library(NMF))

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
comp <- comp[rowSums(comp) > 0, colSums(comp) > 0, drop = FALSE]
meta <- meta[match(rownames(comp), meta$unit_id), , drop = FALSE]

make_posdev <- function(mat, tau = 0.3, scale_factor = NULL) {
  x <- mat
  if (!is.null(scale_factor)) {
    x <- log1p(x * scale_factor)
  }
  z <- scale(x)
  z[!is.finite(z)] <- 0
  out <- pmax(z - tau, 0)
  out[, colSums(out) > 0, drop = FALSE]
}

# Original composition-only block
prop <- sweep(comp, 1, rowSums(comp) + 1e-8, "/")
prop <- prop[, colMeans(prop > 0) >= 0.05, drop = FALSE]
prop_pos <- make_posdev(prop, tau = 0.3, scale_factor = 100)

# Absolute / structural blocks to preserve maturity and density information
pick <- function(name) if (name %in% colnames(comp)) comp[, name] else NULL
sum_available <- function(candidates) {
  keep <- Filter(Negate(is.null), lapply(candidates, pick))
  if (length(keep) == 0) return(NULL)
  Reduce(`+`, keep)
}
abs_feat <- list(
  total_abundance = rowSums(comp),
  lymphoid_total = sum_available(c("B", "Plasma", "CD4_T", "CD8_T", "Tfh-like_CD4", "NK")),
  ilc_total = sum_available(c("ILC1", "ILC2", "ILC3")),
  myeloid_total = sum_available(c("Dendritic", "Macrophage", "Monocyte")),
  b_follicle_axis = sum_available(c("B", "Plasma", "FDC", "Tfh-like_CD4")),
  hev_fdc_axis = sum_available(c("HEV-like_endothelial", "FDC", "B", "Tfh-like_CD4")),
  t_cell_axis = sum_available(c("CD4_T", "CD8_T", "Tfh-like_CD4", "NK")),
  neural_glial_axis = sum_available(c("Glioma", "Glial")),
  vascular_axis = sum_available(c("HEV-like_endothelial", "Vascular")),
  plasma_b_axis = sum_available(c("Plasma", "B"))
)
all_null <- vapply(abs_feat, is.null, logical(1))
if (any(all_null)) {
  warning(sprintf(
    "Dropping engineered axes with no supporting cell types: %s",
    paste(names(abs_feat)[all_null], collapse = ", ")
  ))
  abs_feat <- abs_feat[!all_null]
}
abs_feat <- as.data.frame(lapply(abs_feat, as.numeric), check.names = FALSE)
rownames(abs_feat) <- rownames(comp)
abs_pos <- make_posdev(as.matrix(abs_feat), tau = 0.3, scale_factor = NULL)

x_nmf <- rbind(
  `celltype::` = t(prop_pos),
  `axis::` = t(abs_pos)
)
# The prefixed rownames above are not preserved by rbind; rebuild them explicitly.
rownames(x_nmf) <- c(
  paste0("celltype::", colnames(prop_pos)),
  paste0("axis::", colnames(abs_pos))
)
sample_n <- table(meta$sample)
sample_w <- 1 / as.numeric(sample_n[meta$sample])
names(sample_w) <- meta$unit_id
x_nmf <- sweep(x_nmf, 2, sample_w[colnames(x_nmf)], "*")

cat(sprintf(
  "Augmented NMF input: %d features (%d celltype + %d axis) x %d components\n",
  nrow(x_nmf), ncol(prop_pos), ncol(abs_pos), ncol(x_nmf)
))

set.seed(42)
k <- 5
fit <- nmf(x_nmf, rank = k, method = "brunet", nrun = 50, seed = 42)
W <- basis(fit)
H <- coef(fit)
colnames(W) <- paste0("C", 1:k)
rownames(H) <- paste0("C", 1:k)
colnames(H) <- colnames(x_nmf)

dom <- apply(H, 2, which.max)
dominant <- setNames(paste0("C", dom), colnames(x_nmf))
wt <- as.data.frame(t(H))
wt$unit_id <- rownames(wt)
wt$dominant <- dominant[rownames(wt)]
wt <- wt[, c("unit_id", paste0("C", 1:k), "dominant")]

basis_df <- data.frame(feature = rownames(W), W, row.names = NULL, check.names = FALSE)
write.csv(basis_df, file.path(root, "posdev_aug_tau03_rank5_basis.csv"), row.names = FALSE)
write.csv(wt, file.path(root, "posdev_aug_tau03_rank5_weights.csv"), row.names = FALSE)

summary_df <- do.call(rbind, lapply(seq_len(k), function(i) {
  ord <- order(W[, i], decreasing = TRUE)
  topn <- head(ord, 10)
  data.frame(
    program = paste0("C", i),
    rank = seq_along(topn),
    feature = rownames(W)[topn],
    weight = W[topn, i],
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, file.path(root, "posdev_aug_tau03_rank5_top_features.csv"), row.names = FALSE)

cat("Top features per augmented program:\n")
for (i in seq_len(k)) {
  sub <- summary_df[summary_df$program == paste0("C", i), , drop = FALSE]
  cat(sprintf("%s: %s\n", paste0("C", i), paste(sub$feature[1:5], collapse = ", ")))
}
