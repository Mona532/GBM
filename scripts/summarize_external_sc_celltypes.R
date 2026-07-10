args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript scripts/summarize_external_sc_celltypes.R <out_csv> <file1.rds> [file2.rds ...]")
}

out_csv <- args[1]
files <- args[-1]

rows <- list()

for (f in files) {
  obj <- readRDS(f)
  if (!("Seurat" %in% class(obj))) {
    next
  }

  md <- obj@meta.data
  if (!("Celltype" %in% colnames(md))) {
    next
  }

  tab <- sort(table(md$Celltype), decreasing = TRUE)
  dataset_name <- if ("Dataset" %in% colnames(md)) unique(as.character(md$Dataset)) else basename(f)
  if (length(dataset_name) != 1) {
    dataset_name <- basename(f)
  }

  rows[[length(rows) + 1]] <- data.frame(
    file = basename(f),
    dataset = dataset_name,
    celltype = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)
