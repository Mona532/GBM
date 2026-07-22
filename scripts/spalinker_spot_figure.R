# Generate clean spot-only figures (no tissue background)
library(Seurat)
library(SpaLinker)
library(ggplot2)
library(patchwork)

OUT_DIR  <- "E:/GBM/results/spalinker_tls_v3"
samples_ok <- c("mgh258","zh881inf","zh881t1","zh916bulk","zh1007inf","zh1007nec",
                "zh1019inf","zh1019t1","zh8811a","zh8811b","zh8812")

tumor_region <- c(mgh258="unannotated", zh881inf="infiltrating", zh881t1="T1",
                  zh916bulk="bulk", zh1007inf="infiltrating", zh1007nec="necrotic",
                  zh1019inf="infiltrating", zh1019t1="T1",
                  zh8811a="bulk", zh8811b="bulk", zh8812="bulk")

# ---- Per-sample 4-panel (no tissue image) ----
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next
  se <- readRDS(rds_file)
  region <- tumor_region[[s]]

  pdf(file.path(OUT_DIR, paste0(s, "_spot_only.pdf")), width=14, height=14)

  p1 <- SpatialFeaturePlot(se, features="TLS_score", pt.size.factor=2.0,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title=paste(s, "(", region, ")", "- TLS Score"),
         subtitle=sprintf("max=%.3f  median=%.3f", max(se$TLS_score,na.rm=TRUE), median(se$TLS_score,na.rm=TRUE))) +
    scale_fill_viridis_c(option="inferno")
  p2 <- SpatialFeaturePlot(se, features="B_Plasma", pt.size.factor=2.0,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title="B/Plasma Score") + scale_fill_viridis_c(option="magma")
  p3 <- SpatialFeaturePlot(se, features="T_cell", pt.size.factor=2.0,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title="T Cell Score") + scale_fill_viridis_c(option="magma")
  p4 <- SpatialFeaturePlot(se, features="BT_codist", pt.size.factor=2.0,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title="B-T Co-distribution",
         subtitle=sprintf("max=%.3f", max(se$BT_codist,na.rm=TRUE))) +
    scale_fill_viridis_c(option="plasma")

  print((p1 | p2) / (p3 | p4))
  dev.off()
  cat(sprintf("  %s done\n", s))
}

# ---- Combined overview: TLS score (spot only, all 11 samples) ----
cat("Building overview...\n")
pdf(file.path(OUT_DIR, "GBM_TLS_overview_spot.pdf"), width=22, height=28)
plots <- list()
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next
  se <- readRDS(rds_file)
  region <- tumor_region[[s]]
  plots[[s]] <- SpatialFeaturePlot(se, features="TLS_score", pt.size.factor=1.4,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title=sprintf("%s (%s)", s, region),
         subtitle=sprintf("max=%.3f", max(se$TLS_score,na.rm=TRUE))) +
    scale_fill_viridis_c(option="inferno")
}
print(wrap_plots(plots, ncol=4))
dev.off()

# ---- Combined overview: B-T co-distribution (spot only) ----
pdf(file.path(OUT_DIR, "GBM_BT_codist_overview_spot.pdf"), width=22, height=28)
plots <- list()
for(s in samples_ok) {
  rds_file <- file.path(OUT_DIR, paste0(s, "_seurat.rds"))
  if(!file.exists(rds_file)) next
  se <- readRDS(rds_file)
  region <- tumor_region[[s]]
  plots[[s]] <- SpatialFeaturePlot(se, features="BT_codist", pt.size.factor=1.4,
        image.alpha=0, stroke=0) +
    ggplot2::labs(title=sprintf("%s (%s)", s, region),
         subtitle=sprintf("max=%.3f", max(se$BT_codist,na.rm=TRUE))) +
    scale_fill_viridis_c(option="plasma")
}
print(wrap_plots(plots, ncol=4))
dev.off()

cat("Done!\n")
