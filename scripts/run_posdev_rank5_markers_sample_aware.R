library(edgeR)
library(limma)

root <- "E:/GBM/results"
counts_path <- file.path(root, "tls_pseudobulk_counts_by_component.rds")
weights_path <- file.path(root, "posdev_tau03_rank5_weights.csv")
basis_path <- file.path(root, "posdev_tau03_rank5_basis.csv")

counts <- readRDS(counts_path)
weights <- read.csv(weights_path, check.names = FALSE, stringsAsFactors = FALSE)
basis <- read.csv(basis_path, check.names = FALSE, stringsAsFactors = FALSE)
colnames(basis)[1] <- "cell_type"

weights <- weights[match(colnames(counts), weights$unit_id), , drop = FALSE]
if (any(is.na(weights$unit_id))) {
  stop("Failed to align posdev weights to count matrix columns.")
}

weights$sample <- sub("__c.*$", "", weights$unit_id)
weights$component_id <- sub("^.*__c", "", weights$unit_id)
weights$dominant <- factor(weights$dominant, levels = paste0("C", 1:5))
weights$dominant_weight <- mapply(function(dom, i) {
  as.numeric(weights[[as.character(dom)]][i])
}, weights$dominant, seq_len(nrow(weights)))

sample_program_key <- paste(weights$sample, weights$dominant, sep = "__")
split_idx <- split(seq_len(ncol(counts)), sample_program_key)

agg_counts <- do.call(cbind, lapply(split_idx, function(idx) {
  rowSums(counts[, idx, drop = FALSE])
}))
agg_meta <- do.call(rbind, lapply(names(split_idx), function(key) {
  idx <- split_idx[[key]]
  sample_id <- weights$sample[idx[1]]
  program <- as.character(weights$dominant[idx[1]])
  data.frame(
    agg_id = key,
    sample = sample_id,
    program = program,
    n_components = length(idx),
    mean_dominant_weight = mean(weights$dominant_weight[idx]),
    median_dominant_weight = median(weights$dominant_weight[idx]),
    stringsAsFactors = FALSE
  )
}))
colnames(agg_counts) <- agg_meta$agg_id
agg_meta$program <- factor(agg_meta$program, levels = paste0("C", 1:5))

write.csv(agg_meta, file.path(root, "posdev_tau03_rank5_sample_program_metadata.csv"), row.names = FALSE)

program_summary <- do.call(rbind, lapply(levels(agg_meta$program), function(program) {
  sub <- agg_meta[agg_meta$program == program, , drop = FALSE]
  basis_col <- basis[[program]]
  ord <- order(basis_col, decreasing = TRUE)
  top_basis <- basis$cell_type[ord][1:5]
  data.frame(
    program = program,
    n_sample_program = nrow(sub),
    n_samples = length(unique(sub$sample)),
    mean_n_components = mean(sub$n_components),
    median_n_components = median(sub$n_components),
    mean_dominant_weight = mean(sub$mean_dominant_weight),
    top1 = unname(top_basis[1]),
    top2 = unname(top_basis[2]),
    top3 = unname(top_basis[3]),
    top4 = unname(top_basis[4]),
    top5 = unname(top_basis[5]),
    stringsAsFactors = FALSE
  )
}))
write.csv(program_summary, file.path(root, "posdev_tau03_rank5_sample_aware_program_summary.csv"), row.names = FALSE)

group <- agg_meta$program
keep <- filterByExpr(agg_counts, group = group)
y <- DGEList(counts = agg_counts[keep, , drop = FALSE], group = group)
y <- calcNormFactors(y)
design <- model.matrix(~ sample + group, data = agg_meta)
v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)

coef_names <- colnames(design)
group_coefs <- grep("^group", coef_names, value = TRUE)

make_program_contrast <- function(target) {
  v <- rep(0, length(coef_names))
  names(v) <- coef_names
  others <- setdiff(levels(group), target)
  if (target == "C1") {
    for (other in others) v[paste0("group", other)] <- -1 / length(others)
  } else {
    v[paste0("group", target)] <- 1
    nonbaseline_others <- setdiff(others, "C1")
    for (other in nonbaseline_others) v[paste0("group", other)] <- -1 / length(others)
  }
  v
}

for (program in levels(group)) {
  contrast <- make_program_contrast(program)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  tab <- topTable(fit2, number = Inf, sort.by = "P")
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
  write.csv(tab, file.path(root, paste0("posdev_tau03_rank5_sample_aware_", program, "_markers_vs_rest.csv")), row.names = FALSE)

  top_up <- tab[tab$logFC > 0, c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")]
  write.csv(head(top_up, 30), file.path(root, paste0("posdev_tau03_rank5_sample_aware_", program, "_top30_up.csv")), row.names = FALSE)
}

cat("Saved sample-aware posdev rank5 marker tables\n")
print(program_summary)
