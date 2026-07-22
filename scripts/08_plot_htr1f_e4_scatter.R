suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

root <- "E:/GBM/results"
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
comp <- readRDS(file.path(root, "tls_pseudobulk_c2l_by_component.rds"))
weights <- read.csv(file.path(root, "tls_compnmf_rank5_unit_weights.csv"), check.names = FALSE, stringsAsFactors = FALSE)

weights <- weights[match(colnames(counts), weights$unit_id), , drop = FALSE]
if (any(is.na(weights$unit_id)) || !identical(weights$unit_id, colnames(counts))) {
  stop("Failed to align weights with component counts.")
}

gene <- "HTR1F"
if (!gene %in% rownames(counts)) {
  stop("HTR1F not found in component count matrix.")
}

logcpm <- cpm(DGEList(counts = counts), log = TRUE, prior.count = 2)
e4_units <- weights$unit_id[weights$dominant_ecotype == "E4"]
if (length(e4_units) < 3) {
  stop("Too few E4 components.")
}

plot_df <- do.call(rbind, lapply(c("ILC1", "ILC2", "ILC3"), function(ct) {
  data.frame(
    unit_id = e4_units,
    ilc_type = ct,
    ilc_abundance = as.numeric(comp[e4_units, ct]),
    receptor_expr = as.numeric(logcpm[gene, e4_units]),
    stringsAsFactors = FALSE
  )
}))

stats_df <- do.call(rbind, lapply(split(plot_df, plot_df$ilc_type), function(df) {
  ct <- suppressWarnings(cor.test(df$ilc_abundance, df$receptor_expr, method = "spearman", exact = FALSE))
  data.frame(
    ilc_type = df$ilc_type[1],
    rho = unname(ct$estimate),
    pvalue = ct$p.value,
    label = sprintf("rho = %.2f\nP = %.3g", unname(ct$estimate), ct$p.value),
    stringsAsFactors = FALSE
  )
}))

plot_df$ilc_type <- factor(plot_df$ilc_type, levels = c("ILC1", "ILC2", "ILC3"))
stats_df$ilc_type <- factor(stats_df$ilc_type, levels = c("ILC1", "ILC2", "ILC3"))

p <- ggplot(plot_df, aes(x = ilc_abundance, y = receptor_expr)) +
  geom_point(shape = 21, size = 2.3, stroke = 0.3, colour = "black", fill = "#C97A63", alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, colour = "red4") +
  facet_wrap(~ ilc_type, scales = "free_x", nrow = 1) +
  geom_text(
    data = stats_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.6,
    colour = "black"
  ) +
  labs(
    title = "HTR1F vs ILC abundance in E4 components",
    x = "Component ILC abundance",
    y = "HTR1F expression (logCPM)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(root, "fig_tls_rank5_E4_HTR1F_scatter.jpg"), p, width = 10, height = 3.8, dpi = 300)
ggsave(file.path(root, "fig_tls_rank5_E4_HTR1F_scatter.pdf"), p, width = 10, height = 3.8)
write.csv(plot_df, file.path(root, "tls_rank5_E4_HTR1F_scatter_data.csv"), row.names = FALSE)

cat("Saved HTR1F E4 scatter plot\n")
