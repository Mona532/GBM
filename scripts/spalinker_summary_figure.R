# Generate combined summary figure for all GBM TLS analysis
library(Seurat)
library(SpaLinker)
library(ggplot2)
library(patchwork)

OUT_DIR  <- "E:/GBM/results/spalinker_tls_v3"
samples_ok <- c("mgh258","zh881inf","zh881t1","zh916bulk","zh1007inf","zh1007nec",
                "zh1019inf","zh1019t1","zh8811a","zh8811b","zh8812")

# ---- Per-sample detailed PDF ----
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next

  se <- readRDS(rds_file)

  pdf(file.path(OUT_DIR, paste0(s, "_TLS_detailed.pdf")), width=16, height=12)

  p1 <- SpatialFeaturePlot(se, features="TLS_score", pt.size.factor=1.6) +
    ggplot2::labs(title=paste(s, "- TLS Score"), subtitle=paste("Max:", round(max(se$TLS_score,na.rm=TRUE),3)))
  p2 <- SpatialFeaturePlot(se, features="B_Plasma", pt.size.factor=1.6) +
    ggplot2::labs(title=paste(s, "- B/Plasma Score"))
  p3 <- SpatialFeaturePlot(se, features="T_cell", pt.size.factor=1.6) +
    ggplot2::labs(title=paste(s, "- T Cell Score"))
  p4 <- SpatialFeaturePlot(se, features="BT_codist", pt.size.factor=1.6) +
    ggplot2::labs(title=paste(s, "- B-T Co-distribution"))
  p5 <- SpatialFeaturePlot(se, features="LC50sig", pt.size.factor=1.6) +
    ggplot2::labs(title=paste(s, "- LC.50sig"))
  se$TLS_binary <- factor(ifelse(se$TLS_region=="TLS","TLS","Non-TLS"), levels=c("TLS","Non-TLS"))
  p6 <- SpatialDimPlot(se, group.by="TLS_binary", pt.size.factor=1.6,
        cols=c("TLS"="red","Non-TLS"="grey80")) +
    ggplot2::labs(title=paste(s, "- TLS Regions (0 spots)"))

  print((p1 | p2) / (p3 | p4) / (p5 | p6))
  dev.off()
  cat(sprintf("  %s detailed PDF done\n", s))
}

# ---- Combined overview: TLS score across all samples ----
cat("Building combined overview...\n")
pdf(file.path(OUT_DIR, "GBM_TLS_overview.pdf"), width=20, height=24)
plots <- list()
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next
  se <- readRDS(rds_file)
  plots[[s]] <- SpatialFeaturePlot(se, features="TLS_score", pt.size.factor=1.2) +
    ggplot2::labs(title=s, subtitle=sprintf("TLS max=%.3f", max(se$TLS_score, na.rm=TRUE)))
}
print(wrap_plots(plots, ncol=4))
dev.off()

# ---- Combined B-T co-distribution overview ----
pdf(file.path(OUT_DIR, "GBM_BT_codist_overview.pdf"), width=20, height=24)
plots <- list()
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next
  se <- readRDS(rds_file)
  plots[[s]] <- SpatialFeaturePlot(se, features="BT_codist", pt.size.factor=1.2) +
    ggplot2::labs(title=s, subtitle=sprintf("BT codist max=%.3f", max(se$BT_codist,na.rm=TRUE)))
}
print(wrap_plots(plots, ncol=4))
dev.off()

cat("Done! Files in", OUT_DIR, "\n")
