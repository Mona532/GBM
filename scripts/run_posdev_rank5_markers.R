library(edgeR)

root <- "E:/GBM/results"
counts_path <- file.path(root, "tls_pseudobulk_counts_by_component.rds")
weights_path <- file.path(root, "posdev_tau03_rank5_weights.csv")
basis_path <- file.path(root, "posdev_tau03_rank5_basis.csv")

counts <- readRDS(counts_path)
weights <- read.csv(weights_path, check.names = FALSE, stringsAsFactors = FALSE)
basis <- read.csv(basis_path, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(!is.null(colnames(counts)))
stopifnot("unit_id" %in% names(weights), "dominant" %in% names(weights))

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

basis_rownames <- names(basis)[1]
colnames(basis)[1] <- "cell_type"

program_summary <- do.call(rbind, lapply(levels(weights$dominant), function(program) {
  sub <- weights[weights$dominant == program, , drop = FALSE]
  basis_col <- basis[[program]]
  ord <- order(basis_col, decreasing = TRUE)
  top_basis <- basis$cell_type[ord][1:5]
  data.frame(
    program = program,
    n_units = nrow(sub),
    n_samples = length(unique(sub$sample)),
    mean_dominant_weight = unname(mean(sub$dominant_weight)),
    median_dominant_weight = unname(median(sub$dominant_weight)),
    top1 = unname(top_basis[1]),
    top2 = unname(top_basis[2]),
    top3 = unname(top_basis[3]),
    top4 = unname(top_basis[4]),
    top5 = unname(top_basis[5]),
    stringsAsFactors = FALSE
  )
}))
write.csv(program_summary, file.path(root, "posdev_tau03_rank5_program_summary.csv"), row.names = FALSE)

group <- weights$dominant
keep <- filterByExpr(counts, group = group)
y <- DGEList(counts = counts[keep, , drop = FALSE], group = group)
y <- calcNormFactors(y)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

for (program in levels(group)) {
  others <- setdiff(levels(group), program)
  contrast <- rep(-1 / length(others), length(levels(group)))
  names(contrast) <- levels(group)
  contrast[program] <- 1

  qlf <- glmQLFTest(fit, contrast = contrast)
  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
  write.csv(tab, file.path(root, paste0("posdev_tau03_rank5_", program, "_markers_vs_rest.csv")), row.names = FALSE)

  top_up <- tab[tab$logFC > 0, c("gene", "logFC", "logCPM", "PValue", "FDR")]
  write.csv(head(top_up, 30), file.path(root, paste0("posdev_tau03_rank5_", program, "_top30_up.csv")), row.names = FALSE)
}

cat("Saved posdev rank5 program summaries and marker tables\n")
print(program_summary)
