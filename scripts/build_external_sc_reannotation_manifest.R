args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript scripts/build_external_sc_reannotation_manifest.R <out_csv> <file1.rds> [file2.rds ...]")
}

out_csv <- args[1]
files <- args[-1]

rows <- list()

for (f in files) {
  obj <- readRDS(f)
  if (!("Seurat" %in% class(obj))) {
    rows[[length(rows) + 1]] <- data.frame(
      file = basename(f),
      object_class = paste(class(obj), collapse = "|"),
      n_cells = NA_integer_,
      n_genes = NA_integer_,
      assays = NA_character_,
      default_assay = NA_character_,
      counts_present = NA,
      data_present = NA,
      scale_present = NA,
      reductions = NA_character_,
      meta_columns = NA_character_,
      sample_columns_present = NA_character_,
      reannotate_mode = "unsupported",
      stringsAsFactors = FALSE
    )
    next
  }

  assays <- names(obj@assays)
  default_assay <- DefaultAssay(obj)
  assay_obj <- obj@assays[[default_assay]]
  slots <- slotNames(assay_obj)
  md_cols <- colnames(obj@meta.data)
  sample_cols <- c("Dataset", "Sample", "Patient", "orig.ident", "source")
  sample_cols <- sample_cols[sample_cols %in% md_cols]

  rows[[length(rows) + 1]] <- data.frame(
    file = basename(f),
    object_class = paste(class(obj), collapse = "|"),
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    assays = paste(assays, collapse = "|"),
    default_assay = default_assay,
    counts_present = "counts" %in% slots,
    data_present = "data" %in% slots,
    scale_present = "scale.data" %in% slots,
    reductions = paste(names(obj@reductions), collapse = "|"),
    meta_columns = paste(head(md_cols, 20), collapse = "|"),
    sample_columns_present = paste(sample_cols, collapse = "|"),
    reannotate_mode = "drop_old_Celltype_and_reintegrate",
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)
