args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("usage: Rscript scripts/inspect_external_sc_rds.R <file1.rds> [file2.rds ...]")
}

inspect_seurat <- function(obj) {
  md <- colnames(obj@meta.data)
  cat("NCELL\t", ncol(obj), "\n", sep = "")
  cat("NGENE\t", nrow(obj), "\n", sep = "")
  cat("META_COLS\t", paste(head(md, 30), collapse = "|"), "\n", sep = "")
}

inspect_sce <- function(obj) {
  md <- colnames(SummarizedExperiment::colData(obj))
  cat("NCELL\t", ncol(obj), "\n", sep = "")
  cat("NGENE\t", nrow(obj), "\n", sep = "")
  cat("META_COLS\t", paste(head(md, 30), collapse = "|"), "\n", sep = "")
}

for (f in args) {
  cat("FILE\t", normalizePath(f, winslash = "/"), "\n", sep = "")
  obj <- readRDS(f)
  cls <- paste(class(obj), collapse = "|")
  cat("CLASS\t", cls, "\n", sep = "")

  if ("Seurat" %in% class(obj)) {
    inspect_seurat(obj)
  } else if ("SingleCellExperiment" %in% class(obj)) {
    inspect_sce(obj)
  } else if (is.list(obj)) {
    cat("LIST_NAMES\t", paste(head(names(obj), 30), collapse = "|"), "\n", sep = "")
  } else if (is.matrix(obj) || inherits(obj, "dgCMatrix")) {
    cat("DIMS\t", paste(dim(obj), collapse = "x"), "\n", sep = "")
  } else {
    cat("SUMMARY\t", paste(utils::capture.output(utils::str(obj, max.level = 1)), collapse = " "), "\n", sep = "")
  }

  cat("---\n")
}
