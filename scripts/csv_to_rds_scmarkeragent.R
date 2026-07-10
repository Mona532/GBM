args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript csv_to_rds_scmarkeragent.R <input_csv_gz> <output_rds>")
}

infile <- args[[1]]
outfile <- args[[2]]

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
saveRDS(df, outfile)
cat("saved", outfile, "rows", nrow(df), "cols", ncol(df), "\n")
