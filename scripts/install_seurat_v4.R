options(repos = c(CRAN = "https://cloud.r-project.org"))

lib_dir <- normalizePath("E:/GBM/R/R-4.3.2/library", winslash = "/", mustWork = FALSE)
if (!dir.exists(lib_dir)) {
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
}

.libPaths(c(lib_dir, .libPaths()))

required_bootstrap <- c("remotes")
missing_bootstrap <- setdiff(required_bootstrap, rownames(installed.packages()))
if (length(missing_bootstrap) > 0) {
  install.packages(missing_bootstrap, lib = lib_dir)
}

seurat_related <- c(
  "SeuratObject",
  "Matrix",
  "future",
  "future.apply",
  "ggplot2",
  "patchwork",
  "cowplot",
  "dplyr",
  "data.table",
  "tibble",
  "Rcpp",
  "RcppAnnoy",
  "reticulate",
  "hdf5r",
  "uwot",
  "sctransform",
  "irlba",
  "fitdistrplus",
  "pbapply",
  "progressr",
  "clustree"
)

install.packages(seurat_related, dependencies = TRUE, lib = lib_dir)
remotes::install_version(
  "Seurat",
  version = "4.4.0",
  dependencies = TRUE,
  upgrade = "never",
  lib = lib_dir
)
remotes::install_github(
  "chris-mcginnis-ucsf/DoubletFinder",
  dependencies = TRUE,
  upgrade = "never",
  lib = lib_dir
)

installed <- installed.packages(lib.loc = lib_dir)
targets <- c("Seurat", "SeuratObject", "sctransform", "DoubletFinder", "future", "patchwork")
summary_df <- data.frame(
  package = targets,
  version = vapply(
    targets,
    function(pkg) if (pkg %in% rownames(installed)) installed[pkg, "Version"] else NA_character_,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "E:/GBM/R/install-summary.csv", row.names = FALSE)
print(summary_df)
