options(
  repos = c(CRAN = "https://cran.r-project.org"),
  download.file.method = "libcurl",
  pkgType = "both"
)

lib_dir <- normalizePath("E:/GBM/R/R-4.3.2/library", winslash = "/", mustWork = FALSE)
if (!dir.exists(lib_dir)) {
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
}

.libPaths(c(lib_dir, .libPaths()))

bootstrap <- c("remotes")
missing_bootstrap <- setdiff(bootstrap, rownames(installed.packages()))
if (length(missing_bootstrap) > 0) {
  install.packages(missing_bootstrap, lib = lib_dir)
}

# Keep this lean: install only hard runtime dependencies for Seurat v4.
core_pkgs <- c(
  "Matrix",
  "future",
  "future.apply",
  "ggplot2",
  "patchwork",
  "cowplot",
  "dplyr",
  "Rcpp",
  "RcppAnnoy",
  "reticulate",
  "hdf5r",
  "uwot",
  "sctransform",
  "irlba"
)

missing_core <- setdiff(core_pkgs, rownames(installed.packages()))
if (length(missing_core) > 0) {
  install.packages(
    missing_core,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    lib = lib_dir
  )
}

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  remotes::install_version(
    "SeuratObject",
    version = "4.1.4",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    lib = lib_dir
  )
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_version(
    "Seurat",
    version = "4.4.0",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    lib = lib_dir
  )
}

installed <- installed.packages(lib.loc = lib_dir)
targets <- c("Seurat", "SeuratObject", "sctransform", "ggplot2", "future", "patchwork", "uwot")
summary_df <- data.frame(
  package = targets,
  version = vapply(
    targets,
    function(pkg) if (pkg %in% rownames(installed)) installed[pkg, "Version"] else NA_character_,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "E:/GBM/R/install-summary-minimal.csv", row.names = FALSE)
print(summary_df)
