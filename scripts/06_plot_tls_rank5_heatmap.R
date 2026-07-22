suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "E:/GBM/results"
basis <- read.csv(file.path(root, "tls_compnmf_rank5_basis.csv"), check.names = FALSE)
basis_raw <- read.csv(file.path(root, "tls_compnmf_rank5_basis_raw.csv"), check.names = FALSE)
summary_df <- read.csv(file.path(root, "tls_compnmf_rank5_ecotype_summary.csv"), check.names = FALSE)
weight_df <- read.csv(file.path(root, "tls_compnmf_rank5_unit_weights.csv"), check.names = FALSE, stringsAsFactors = FALSE)
comp <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))

rownames(basis) <- basis$cell_type
basis <- as.matrix(basis[, grep("^E[0-9]+$", colnames(basis)), drop = FALSE])
storage.mode(basis) <- "numeric"
rownames(basis_raw) <- basis_raw$cell_type
basis_raw <- as.matrix(basis_raw[, grep("^[0-9]+$", colnames(basis_raw)), drop = FALSE])
storage.mode(basis_raw) <- "numeric"
eco_order <- colnames(basis)
colnames(basis_raw) <- eco_order

# NMF heatmap: row-normalized view to emphasize subtype specificity
basis_prop <- sweep(basis_raw, 2, colSums(basis_raw) + 1e-8, "/")
basis_z <- t(scale(t(basis_prop)))
basis_z[!is.finite(basis_z)] <- 0

counts <- setNames(summary_df$n_units, summary_df$ecotype)

stopifnot(all(weight_df$unit_id %in% rownames(comp)))
comp <- comp[weight_df$unit_id, , drop = FALSE]
weight_df$dominant_ecotype <- factor(weight_df$dominant_ecotype, levels = eco_order)
raw_mean <- sapply(eco_order, function(eco) {
  idx <- weight_df$dominant_ecotype == eco
  colMeans(comp[idx, , drop = FALSE], na.rm = TRUE)
})
raw_mean <- t(raw_mean)
raw_mean <- t(raw_mean)

eco_labels <- setNames(sprintf("%s\n(n=%d)", eco_order, counts[eco_order]), eco_order)
cell_order <- rownames(basis_prop)[order(apply(basis_prop, 1, which.max), -apply(basis_prop, 1, max))]

make_df <- function(mat, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("cell_type", "ecotype", value_name)
  df$cell_type <- factor(df$cell_type, levels = rev(cell_order))
  df$ecotype <- factor(df$ecotype, levels = eco_order, ordered = TRUE)
  df
}

df_raw <- make_df(raw_mean[cell_order, eco_order, drop = FALSE], "value")
df_z <- make_df(basis_z[cell_order, eco_order, drop = FALSE], "value")

p1 <- ggplot(df_raw, aes(x = ecotype, y = cell_type, fill = value)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradient(low = "#f7fbff", high = "#b2182b", name = "Mean abundance") +
  scale_x_discrete(labels = eco_labels) +
  labs(x = NULL, y = NULL, title = "Raw abundance heatmap") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

p2 <- ggplot(df_z, aes(x = ecotype, y = cell_type, fill = value)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "#f7f7f7",
    high = "#b2182b",
    midpoint = 0,
    name = "Score"
  ) +
  scale_x_discrete(labels = eco_labels) +
  labs(x = NULL, y = NULL, title = "NMF heatmap") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

jpg_path1 <- file.path(root, "fig_tls_rank5_abundance_heatmap.jpg")
pdf_path1 <- file.path(root, "fig_tls_rank5_abundance_heatmap.pdf")
jpg_path2 <- file.path(root, "fig_tls_rank5_nmf_heatmap.jpg")
pdf_path2 <- file.path(root, "fig_tls_rank5_nmf_heatmap.pdf")

jpeg(jpg_path1, width = 1200, height = 1200, res = 220)
print(p1)
dev.off()

pdf(pdf_path1, width = 5.8, height = 5.4)
print(p1)
dev.off()

jpeg(jpg_path2, width = 1200, height = 1200, res = 220)
print(p2)
dev.off()

pdf(pdf_path2, width = 5.8, height = 5.4)
print(p2)
dev.off()

cat("Saved fig_tls_rank5_abundance_heatmap and fig_tls_rank5_nmf_heatmap\n")
