# ============================================================
# 02_deg_analysis.R — TLS program differential expression
# ============================================================
# Input:  tls_pseudobulk_counts_by_component.rds (cleaned pseudobulk counts)
#         posdev_tau03_rank5_weights.csv (program assignments)
# Output: posdev_all_DEGs.csv (per-program DEGs with logFC, adj.P.Val)
#         posdev_GO_relaxed.csv (GO-BP enrichment for programs with ≥20 DEGs)
# Method: voom → duplicateCorrelation(block=sample) → lmFit(~0+prog+batch)
#         → contrasts.fit(one-vs-rest) → eBayes(robust=TRUE)
#         DEG threshold: FDR < 0.1 & logFC > 0.5
#         GO via clusterProfiler, Hallmark via msigdbr
#
# Key design decisions:
# - batch as covariate in design matrix (not removeBatchEffect pre-modeling)
# - duplicateCorrelation handles pseudoreplication (multiple components/sample)
# - eBayes AFTER contrasts.fit (standard limma order)
# - filterByExpr for gene-level filtering (preserves program-specific markers)
# - Hallmark uses gene_symbol (not ENTREZID) for consistent ID matching
# ============================================================

suppressMessages({
  library(edgeR)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

root   <- "E:/GBM/results"
counts <- readRDS(file.path(root, "tls_pseudobulk_counts_by_component.rds"))
wt     <- read.csv(file.path(root, "posdev_tau03_rank5_weights.csv"),
                   check.names = FALSE)

# program identity (transcription-driven names, matched to DEG evidence)
idmap <- c(
  C1 = "Inflammatory",
  C2 = "Myeloid",
  C3 = "Neuronal",
  C4 = "Oligodendrocyte",
  C5 = "HypoxicVascular"
)

# ---- prepare data ----
counts <- as.matrix(counts[, wt$unit_id, drop = FALSE])
sample <- sub("__c[0-9]+$", "", colnames(counts))
batch  <- ifelse(grepl("^AT15", sample), "AT15",
          ifelse(grepl("^AT",   sample), "AT",
          ifelse(grepl("^dryad", sample), "dryad", "other_small")))
prog   <- factor(as.character(wt$dominant), levels = c("C1", "C2", "C3", "C4", "C5"))

# ---- design matrix: program + batch as covariates ----
design <- model.matrix(~ 0 + prog + factor(batch))
colnames(design) <- make.names(colnames(design))

# ---- voom pipeline with pseudoreplication correction ----
y    <- DGEList(counts = counts)
keep <- filterByExpr(y, design = design)          # gene-level expression filter
y    <- y[keep, , keep.lib.sizes = FALSE]
y    <- calcNormFactors(y)
cat(sprintf("genes after filterByExpr: %d / %d\n", nrow(y), nrow(counts)))

v  <- voom(y, design)
cf <- duplicateCorrelation(v, design, block = sample)
cat(sprintf("intra-sample correlation: %.3f\n", cf$consensus.correlation))
v  <- voom(y, design, block = sample, correlation = cf$consensus.correlation)

fit <- lmFit(v, design, block = sample, correlation = cf$consensus.correlation)
fit$genes <- data.frame(SYMBOL = rownames(v), stringsAsFactors = FALSE)

# ---- one-vs-rest contrasts for all 5 programs ----
progs <- names(idmap)
pcols <- paste0("prog", c("C1", "C2", "C3", "C4", "C5"))
CM    <- matrix(0, ncol(design), 5, dimnames = list(colnames(design), progs))
for (i in 1:5) {
  CM[pcols[i], i] <- 1
  CM[setdiff(pcols, pcols[i]), i] <- -1/4     # vs mean of other 4
}
fit2 <- contrasts.fit(fit, CM)
fit2 <- eBayes(fit2, robust = TRUE)

# ---- extract DEGs and run enrichment ----
ks     <- keys(org.Hs.eg.db, keytype = "SYMBOL")
uv_sym <- rownames(v)[rownames(v) %in% ks]
uv     <- bitr(uv_sym, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db, drop = FALSE)$ENTREZID
hm_all <- msigdbr(species = "Homo sapiens", collection = "H")   # Hallmark gene sets

all_deg  <- list()
go_all   <- list()
hm_all_res <- list()

for (i in 1:5) {
  p  <- progs[i]
  tt <- topTable(fit2, coef = i, number = Inf)

  # DEGs at relaxed threshold
  sig <- tt[tt$adj.P.Val < 0.1 & tt$logFC > 0.5,
            c("SYMBOL", "logFC", "adj.P.Val", "AveExpr")]
  if (nrow(sig) > 0) { sig$program <- p; all_deg[[p]] <- sig }
  cat(sprintf("%-16s: %d DEGs (FDR<0.1, logFC>0.5)\n", p, nrow(sig)))

  # enrichment only if ≥20 DEGs
  if (nrow(sig) < 20) next
  mapped <- sig$SYMBOL[sig$SYMBOL %in% ks]
  m <- bitr(mapped, fromType = "SYMBOL", toType = "ENTREZID",
            OrgDb = org.Hs.eg.db, drop = FALSE)
  m <- m[!is.na(m$ENTREZID), ]
  if (nrow(m) < 10) next

  # GO-BP
  go <- tryCatch(
    simplify(enrichGO(m$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
              ont = "BP", universe = uv, pAdjustMethod = "BH",
              qvalueCutoff = 0.1, minGSSize = 15, maxGSSize = 500),
             cutoff = 0.7),
    error = function(e) NULL)
  if (!is.null(go) && nrow(go) > 0) {
    d <- as.data.frame(go); d$program <- p; go_all[[p]] <- d
    cat(sprintf("  GO-BP: %d terms. top3: %s\n", nrow(d),
                paste(head(d$Description, 3), collapse = " | ")))
  }

  # Hallmark (gene_symbol matching to avoid ENTREZID/symbol mismatch)
  hm_sub <- hm_all[hm_all$gene_symbol %in% mapped, c("gs_name", "gene_symbol")]
  hr <- tryCatch(
    enricher(mapped, TERM2GENE = hm_sub, universe = uv_sym,
             pAdjustMethod = "BH", qvalueCutoff = 0.25),
    error = function(e) NULL)
  if (!is.null(hr) && nrow(hr) > 0) {
    d <- as.data.frame(hr); d$program <- p; hm_all_res[[p]] <- d
    cat(sprintf("  Hallmark: %d terms. top3: %s\n", nrow(d),
                paste(head(d$ID, 3), collapse = " | ")))
  }
}

# ---- save ----
deg_out <- do.call(rbind, all_deg)
deg_out <- deg_out[order(deg_out$program, -deg_out$logFC), ]
write.csv(deg_out, file.path(root, "posdev_all_DEGs.csv"), row.names = FALSE)

if (length(go_all))
  write.csv(do.call(rbind, go_all),
            file.path(root, "posdev_GO_relaxed.csv"), row.names = FALSE)
if (length(hm_all_res))
  write.csv(do.call(rbind, hm_all_res),
            file.path(root, "posdev_hallmark_relaxed.csv"), row.names = FALSE)

cat(sprintf("\nSaved: %d DEGs total. Done.\n", nrow(deg_out)))
