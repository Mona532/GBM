library(rhdf5)

fp <- "E:/GBM/spatial_data_visium/spatial_data_visium/anndata/AT3-BRA5-FO-1_1.h5ad"

h5ls(fp, recursive=FALSE)
cat("\n--- var group ---\n")
h5ls(fp, "var")

cat("\n--- var/feature_type ---\n")
ft <- h5read(fp, "var/feature_type")
print(table(ft))

# Get index for Cell state abundances
cell_idx <- which(ft == "Cell state abundances")
cat(sprintf("Cell state abundances: %d features\n", length(cell_idx)))

# Get var names
vn <- h5read(fp, "var/_index")
cell_names <- vn[cell_idx]
cat("Cell state names:\n")
for(nm in cell_names) cat(sprintf("  %s\n", nm))

H5close()
