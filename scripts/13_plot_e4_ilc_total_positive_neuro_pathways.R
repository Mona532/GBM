suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(grid)
})

theme_set(
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

root <- "E:/GBM/results"
rank_file <- file.path(root, "e4_spot_gene_corr_ILC_total.csv")
term_file <- file.path(root, "e4_spot_gene_corr_ILC_total_prerank_combined_neuro_hits.csv")
reactome_gmt <- file.path(root, "pathway_library_cache", "Reactome_2022.gmt")
gobp_gmt <- file.path(root, "pathway_library_cache", "GO_Biological_Process_2023.gmt")

out_term_csv <- file.path(root, "e4_spot_gene_corr_ILC_total_positive_neuro_top.csv")
out_source_csv <- file.path(root, "e4_ilc_total_positive_neuro_gsea_source_data.csv")
out_jpg <- file.path(root, "fig_e4_ilc_total_positive_neuro_pathways.jpg")
out_pdf <- file.path(root, "fig_e4_ilc_total_positive_neuro_pathways.pdf")

bad_label <- "^(RPL|RPS|MT-|MTRNR|LINC|MALAT1$|NEAT1$)"

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parsed, function(x) unique(x[-c(1, 2)]))
  names(sets) <- vapply(parsed, `[`, character(1), 1)
  sets
}

clean_term_label <- function(x) {
  x <- gsub(" R-HSA-[0-9]+$", "", x)
  x <- gsub(" \\(GO:[0-9]+\\)$", "", x)
  x <- gsub("^Positive Regulation Of ", "Positive regulation of ", x)
  x <- gsub("^Long-Term Synaptic Potentiation$", "Long-term synaptic potentiation", x)
  x <- gsub("^Synaptic Vesicle Exocytosis$", "Synaptic vesicle exocytosis", x)
  x <- gsub("^Neurotransmitter Secretion$", "Neurotransmitter secretion", x)
  x <- gsub("^Regulation Of Neuronal Synaptic Plasticity$", "Regulation of neuronal synaptic plasticity", x)
  x <- gsub("^Chemical Synaptic Transmission$", "Chemical synaptic transmission", x)
  x <- gsub("^Neurotransmitter Release Cycle$", "Neurotransmitter release cycle", x)
  x <- gsub("^Dopamine Neurotransmitter Release Cycle$", "Dopamine release cycle", x)
  x <- gsub("^Serotonin Neurotransmitter Release Cycle$", "Serotonin release cycle", x)
  x <- gsub("^GABA Synthesis, Release, Reuptake And Degradation$", "GABA synthesis/release/reuptake", x)
  x <- gsub(" Of ", " of ", x, fixed = TRUE)
  x <- gsub(" And ", " and ", x, fixed = TRUE)
  x
}

compute_gsea_curve <- function(ranked_df, gene_set, exponent = 1) {
  genes <- ranked_df$gene
  stats <- ranked_df$rho
  hits <- genes %in% gene_set
  nh <- sum(hits)
  n <- length(genes)
  if (nh < 5 || nh >= n) {
    return(NULL)
  }

  weights <- abs(stats[hits]) ^ exponent
  norm_hit <- sum(weights)
  if (!is.finite(norm_hit) || norm_hit <= 0) {
    return(NULL)
  }

  hit_step <- rep(0, n)
  hit_step[hits] <- (abs(stats[hits]) ^ exponent) / norm_hit
  miss_step <- rep(1 / (n - nh), n)
  miss_step[hits] <- 0
  running <- cumsum(hit_step - miss_step)
  es_idx <- which.max(abs(running))[1]

  data.frame(
    rank = seq_len(n),
    gene = genes,
    rho = stats,
    hit = hits,
    running_score = running,
    es_rank = es_idx,
    stringsAsFactors = FALSE
  )
}

rank_df <- read.csv(rank_file, stringsAsFactors = FALSE)
rank_df <- rank_df[!grepl(bad_label, rank_df$gene, ignore.case = TRUE), c("gene", "rho", "pvalue", "fdr")]
rank_df <- rank_df[complete.cases(rank_df[, c("gene", "rho")]), , drop = FALSE]
rank_df <- rank_df[order(-rank_df$rho, rank_df$pvalue), , drop = FALSE]
rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]

term_df <- read.csv(term_file, check.names = FALSE, stringsAsFactors = FALSE)
term_df <- term_df[term_df$direction == "positive", , drop = FALSE]
term_df$`FDR q-val` <- as.numeric(term_df$`FDR q-val`)
term_df$NES <- as.numeric(term_df$NES)
term_df <- term_df[order(term_df$`FDR q-val`, -term_df$NES), , drop = FALSE]

top_reactome <- head(term_df[term_df$Gene_set == "Reactome_2022", , drop = FALSE], 3)
top_gobp <- head(term_df[term_df$Gene_set == "GO_Biological_Process_2023", , drop = FALSE], 3)
selected_terms <- rbind(top_reactome, top_gobp)
if (nrow(selected_terms) == 0) {
  stop("No positive neuro-related pathways available for plotting.")
}

write.csv(selected_terms, out_term_csv, row.names = FALSE)

gmt_list <- c(read_gmt(reactome_gmt), read_gmt(gobp_gmt))
palette_vals <- c(
  "#24496B", "#4C7DAA", "#7FB0D1",
  "#8E4D2B", "#C87746", "#E4B07A"
)

curve_list <- list()
source_list <- list()
for (i in seq_len(nrow(selected_terms))) {
  term <- selected_terms$Term[i]
  gene_set <- gmt_list[[term]]
  if (is.null(gene_set)) {
    next
  }
  curve <- compute_gsea_curve(rank_df, gene_set)
  if (is.null(curve)) {
    next
  }
  curve$term <- term
  curve$term_label <- clean_term_label(term)
  curve$panel <- ifelse(selected_terms$Gene_set[i] == "Reactome_2022", "Reactome", "GO Biological Process")
  curve$NES <- selected_terms$NES[i]
  curve$FDR <- selected_terms$`FDR q-val`[i]
  curve$color <- palette_vals[i]
  curve_list[[length(curve_list) + 1]] <- curve

  source_list[[length(source_list) + 1]] <- data.frame(
    term = term,
    term_label = clean_term_label(term),
    panel = curve$panel[1],
    rank = curve$rank,
    gene = curve$gene,
    rho = curve$rho,
    hit = curve$hit,
    running_score = curve$running_score,
    NES = curve$NES[1],
    FDR = curve$FDR[1],
    stringsAsFactors = FALSE
  )
}

if (length(curve_list) == 0) {
  stop("Failed to reconstruct any GSEA running score curves.")
}

curve_df <- do.call(rbind, curve_list)
source_df <- do.call(rbind, source_list)
write.csv(source_df, out_source_csv, row.names = FALSE)

term_levels <- unique(curve_df$term_label)
curve_df$term_label <- factor(curve_df$term_label, levels = term_levels)

tick_df <- unique(curve_df[curve_df$hit, c("rank", "term", "term_label", "panel", "color")])
tick_df <- tick_df[order(tick_df$term_label, tick_df$rank), , drop = FALSE]
tick_df$ymin <- as.numeric(tick_df$term_label) - 0.35
tick_df$ymax <- as.numeric(tick_df$term_label) + 0.35

summary_df <- unique(curve_df[, c("term_label", "panel", "NES", "FDR", "color")])
summary_df <- summary_df[order(summary_df$panel, summary_df$term_label), , drop = FALSE]
summary_df$pval <- selected_terms$`NOM p-val`[match(as.character(summary_df$term_label), clean_term_label(selected_terms$Term))]
summary_df$label <- sprintf(
  "%s  NES=%.2f  P=%.3g  FDR=%.3g",
  summary_df$term_label,
  summary_df$NES,
  summary_df$pval,
  summary_df$FDR
)

top_panel <- ggplot(curve_df, aes(rank, running_score, colour = term_label)) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.3, colour = "#6F6F6F") +
  geom_line(linewidth = 0.8, alpha = 0.95) +
  scale_colour_manual(
    values = setNames(summary_df$color, summary_df$term_label)
  ) +
  labs(
    title = "E4 ILC-total positively associated neuro pathways",
    x = NULL,
    y = "Running enrichment score",
    colour = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

mid_panel <- ggplot(tick_df, aes(rank, colour = term_label)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), linewidth = 0.35, alpha = 0.9) +
  scale_colour_manual(values = setNames(summary_df$color, summary_df$term_label), guide = "none") +
  scale_y_continuous(
    breaks = seq_along(term_levels),
    labels = term_levels,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(0, 5.5, 0, 5.5)
  )

rho_df <- unique(curve_df[, c("rank", "rho")])
bottom_panel <- ggplot(rho_df, aes(rank, rho)) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.3, colour = "#6F6F6F") +
  geom_segment(aes(xend = rank, yend = 0), linewidth = 0.12, colour = "#B9B9B9") +
  geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.8, colour = "#222222") +
  labs(x = "Rank in ordered gene list", y = "Spearman rho") +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

legend_df <- summary_df
legend_df$idx <- seq_len(nrow(legend_df))
ncol_legend <- 2
nrow_legend <- ceiling(nrow(legend_df) / ncol_legend)
legend_df$col <- (legend_df$idx - 1) %/% nrow_legend
legend_df$row <- (legend_df$idx - 1) %% nrow_legend
legend_df$x0 <- 0.02 + legend_df$col * 0.50
legend_df$x1 <- legend_df$x0 + 0.07
legend_df$x_text <- legend_df$x1 + 0.02
legend_df$y <- 0.92 - legend_df$row * 0.28

legend_panel <- ggplot(legend_df) +
  geom_segment(aes(x = x0, xend = x1, y = y, yend = y, colour = term_label), linewidth = 1.1) +
  geom_text(aes(x = x_text, y = y, label = label), hjust = 0, vjust = 0.5, size = 3.0) +
  scale_colour_manual(values = setNames(summary_df$color, summary_df$term_label), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

final_plot <- plot_grid(
  top_panel,
  mid_panel,
  bottom_panel,
  legend_panel,
  ncol = 1,
  align = "v",
  rel_heights = c(2.2, 1.1, 1.2, 0.55)
)

ggsave(out_jpg, final_plot, width = 12.5, height = 8.6, dpi = 320, bg = "white")
ggsave(out_pdf, final_plot, width = 12.5, height = 8.6, bg = "white", device = cairo_pdf)

print(summary_df[order(summary_df$panel, -summary_df$NES), c("panel", "term_label", "NES", "FDR")])
