options(width = 220)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(RColorBrewer)
  library(scales)
})

theme_nature_umap <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.size = unit(3.5, "mm"),
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2),
      panel.grid = element_blank()
    )
}

theme_nature_dotplot <- function(base_size = 6.2, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.6),
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      panel.grid = element_blank()
    )
}

save_pub_r <- function(plot_obj, filename, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = width_in, height = height_in)
  print(plot_obj)
  dev.off()

  grDevices::cairo_pdf(
    paste0(filename, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot_obj)
  dev.off()

  ragg::agg_tiff(
    paste0(filename, ".tiff"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  print(plot_obj)
  dev.off()
}

make_palette <- function(labels) {
  n <- length(labels)
  if (n <= 12) {
    vals <- brewer.pal(max(3, n), "Set3")[seq_len(n)]
  } else {
    vals <- grDevices::hcl(
      h = seq(15, 375, length.out = n + 1)[seq_len(n)],
      l = 68,
      c = 80
    )
  }
  names(vals) <- labels
  vals
}

pick_assay <- function(obj) {
  for (assay_name in c("SCT", "RNA", "integrated")) {
    if (assay_name %in% Assays(obj)) {
      DefaultAssay(obj) <- assay_name
      return(obj)
    }
  }
  obj
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

build_umap_plot <- function(obj, group_var, title_text, palette_values, ncol_legend = 1) {
  p <- DimPlot(
    object = obj,
    reduction = "umap",
    group.by = group_var,
    cols = palette_values,
    raster = TRUE,
    shuffle = TRUE,
    pt.size = 0.15
  ) +
    labs(title = title_text, x = "UMAP1", y = "UMAP2", colour = NULL) +
    theme_nature_umap() +
    guides(
      colour = guide_legend(
        override.aes = list(size = 2.5, alpha = 1),
        ncol = ncol_legend,
        byrow = TRUE
      )
    )
  p
}

build_feature_plot <- function(obj, features, title_text) {
  feature_list <- FeaturePlot(
    object = obj,
    reduction = "umap",
    features = features,
    raster = TRUE,
    combine = FALSE,
    order = TRUE,
    pt.size = 0.15,
    min.cutoff = "q05",
    max.cutoff = "q95",
    cols = c("#E8E8E8", "#D24B40")
  )

  feature_list <- lapply(feature_list, function(p) {
    p +
      theme_nature_umap() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 6.4, face = "bold")
      )
  })

  wrap_plots(feature_list, ncol = 4) +
    plot_annotation(title = title_text) &
    theme(plot.title = element_text(size = 7.2, face = "bold"))
}

build_dotplot <- function(obj, features, group_var, title_text, celltype_order) {
  obj[[group_var]] <- factor(obj[[group_var, drop = TRUE]], levels = celltype_order)
  p <- DotPlot(
    object = obj,
    features = features,
    group.by = group_var,
    cols = c("#E5E5E5", "#3182BD"),
    dot.scale = 5.5
  ) +
    scale_colour_gradientn(colours = c("#E5E5E5", "#9ECAE1", "#3182BD", "#08519C")) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL,
      colour = "Average expression",
      size = "Pct. expressing"
    ) +
    theme_nature_dotplot()
  p
}

dataset_specs <- list(
  fetal_liver = list(
    path = "E:/blood_cell_development/geo_data/geo_bundle/fetal_4t_pc/liver/fetal_fetal_liver_sct_cd45.rds",
    label = "Fetal liver CD45+",
    feature_markers = c("CD3D", "MS4A1", "NKG7", "LYZ", "S100A8", "MPO", "HBB", "PPBP", "KIT", "MKI67", "CD74", "JCHAIN")
  ),
  fetal_bm = list(
    path = "E:/blood_cell_development/geo_data/geo_bundle/fetal_4t_pc/bone_marrow/fetal_bm_sct_cd45.rds",
    label = "Fetal bone marrow CD45+",
    feature_markers = c("CD3D", "MS4A1", "NKG7", "LYZ", "S100A8", "MPO", "HBB", "PPBP", "KIT", "MKI67", "CD74", "JCHAIN")
  )
)

output_root <- file.path("results", "fetal_cd45_nature")
ensure_dir(output_root)

for (spec_name in names(dataset_specs)) {
  spec <- dataset_specs[[spec_name]]
  out_dir <- file.path(output_root, sanitize_name(spec_name))
  ensure_dir(out_dir)

  message("Reading ", spec$path)
  obj <- readRDS(spec$path)
  obj <- pick_assay(obj)

  if (!"umap" %in% Reductions(obj)) {
    stop("UMAP reduction not found for ", spec_name)
  }
  if (!all(c("dataset", "celltype") %in% colnames(obj@meta.data))) {
    stop("Expected metadata columns dataset and celltype for ", spec_name)
  }

  meta_df <- obj@meta.data
  meta_df$cell_id <- rownames(meta_df)
  write.csv(meta_df, file.path(out_dir, paste0(spec_name, "_metadata.csv")), row.names = FALSE)

  dataset_levels <- names(sort(table(obj$dataset), decreasing = TRUE))
  celltype_levels <- names(sort(table(obj$celltype), decreasing = TRUE))

  dataset_palette <- make_palette(dataset_levels)
  celltype_palette <- make_palette(celltype_levels)

  write.csv(
    data.frame(group = names(dataset_palette), colour = unname(dataset_palette)),
    file.path(out_dir, paste0(spec_name, "_dataset_palette.csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(group = names(celltype_palette), colour = unname(celltype_palette)),
    file.path(out_dir, paste0(spec_name, "_celltype_palette.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(sort(table(obj$dataset), decreasing = TRUE)),
    file.path(out_dir, paste0(spec_name, "_dataset_counts.csv")),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(sort(table(obj$celltype), decreasing = TRUE)),
    file.path(out_dir, paste0(spec_name, "_celltype_counts.csv")),
    row.names = FALSE
  )

  p_dataset <- build_umap_plot(
    obj = obj,
    group_var = "dataset",
    title_text = paste0(spec$label, " by dataset"),
    palette_values = dataset_palette,
    ncol_legend = ifelse(length(dataset_levels) > 8, 2, 1)
  )
  save_pub_r(
    p_dataset,
    file.path(out_dir, paste0(spec_name, "_dataset_umap")),
    width_mm = 183,
    height_mm = 135
  )

  p_celltype <- build_umap_plot(
    obj = obj,
    group_var = "celltype",
    title_text = paste0(spec$label, " by cell type"),
    palette_values = celltype_palette,
    ncol_legend = ifelse(length(celltype_levels) > 45, 3, ifelse(length(celltype_levels) > 20, 2, 1))
  )
  save_pub_r(
    p_celltype,
    file.path(out_dir, paste0(spec_name, "_celltype_umap")),
    width_mm = ifelse(length(celltype_levels) > 45, 260, 220),
    height_mm = ifelse(length(celltype_levels) > 80, 190, 160)
  )

  feature_markers <- spec$feature_markers[spec$feature_markers %in% rownames(obj)]
  if (length(feature_markers) == 0) {
    stop("No requested feature markers are present in ", spec_name)
  }
  write.csv(
    data.frame(marker = feature_markers),
    file.path(out_dir, paste0(spec_name, "_feature_markers.csv")),
    row.names = FALSE
  )

  p_feature <- build_feature_plot(
    obj = obj,
    features = feature_markers,
    title_text = paste0(spec$label, " classical marker features")
  )
  save_pub_r(
    p_feature,
    file.path(out_dir, paste0(spec_name, "_featureplot_markers")),
    width_mm = 220,
    height_mm = 180
  )

  dot_markers <- feature_markers
  dot_df <- FetchData(obj, vars = c("celltype", dot_markers))
  write.csv(
    dot_df,
    file.path(out_dir, paste0(spec_name, "_dotplot_source_data.csv")),
    row.names = FALSE
  )

  p_dot <- build_dotplot(
    obj = obj,
    features = dot_markers,
    group_var = "celltype",
    title_text = paste0(spec$label, " classical marker dot plot"),
    celltype_order = rev(celltype_levels)
  )
  save_pub_r(
    p_dot,
    file.path(out_dir, paste0(spec_name, "_dotplot_markers")),
    width_mm = 180,
    height_mm = max(120, 22 + 2.5 * length(celltype_levels))
  )
}

message("Done. Outputs written to: ", normalizePath(output_root, winslash = "/"))
