# ============================================================
# 03_deg_heatmap.R — DEG heatmap expression matrix
# ============================================================
# Input:  posdev_all_DEGs.csv (from 02_deg_analysis.R)
#         tls_pseudobulk_counts_by_component.rds
#         posdev_tau03_rank5_weights.csv
# Output: posdev_DEG_heatmap_data.csv (gene × program mean logCPM matrix)
# ============================================================
# Note: only programs with ≥20 significant DEGs are included.
#       Expression values are log2(CPM + 1), program-level means.
#       The DEG filter (FDR<0.1, logFC>0.5) was applied in 02_deg_analysis.R.

suppressMessages(library(edgeR))

root <- "E:/GBM/results"
deg  <- read.csv(file.path(root, "posdev_all_DEGs.csv"))
cnt  <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
wt   <- read.csv(file.path(root, "posdev_tau03_rank5_weights.csv"),
                 check.names = FALSE)

# program mapping: C-code → display name
idmap <- c(C1 = "ILC-NK", C2 = "Myeloid",  C3 = "TLS-structural",
           C4 = "Lymphoid-Glial", C5 = "T-cell")
lab <- setNames(idmap[as.character(wt$dominant)], wt$unit_id)

# all 5 programs, top 10 DEGs each by logFC
progs <- c(C1 = "ILC-NK", C2 = "Myeloid", C3 = "TLS-structural",
           C4 = "Lymphoid-Glial", C5 = "T-cell")
topN <- 10
genes <- unlist(lapply(names(progs), function(p) {
  head(deg$SYMBOL[deg$program == p], topN)
}))
genes <- intersect(genes, rownames(cnt))
prog_vec <- rep(unname(progs), each = topN)
prog_vec <- prog_vec[seq_along(genes)]

# per-program mean logCPM (subset genes only, pre-computed lib sizes)
lib_sizes <- colSums(cnt[, wt$unit_id, drop = FALSE])
y <- DGEList(as.matrix(cnt[genes, wt$unit_id, drop = FALSE]),
             lib.size = lib_sizes)
logcpm <- cpm(calcNormFactors(y), log = TRUE, prior.count = 1)

M <- sapply(progs, function(p) rowMeans(logcpm[genes, lab == p, drop = FALSE]))
colnames(M) <- unname(progs); rownames(M) <- genes

out <- data.frame(gene = genes, program = prog_vec,
                  round(M, 4), check.names = FALSE)
write.csv(out, file.path(root, "posdev_DEG_heatmap_data.csv"),
          row.names = FALSE)

cat(sprintf("%d genes × %d programs, range [%.1f, %.1f], anyNA = %s\n",
            nrow(M), ncol(M), min(M), max(M), any(is.na(M))))
