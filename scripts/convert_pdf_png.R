library(magick)

BASE <- "E:/GBM/results/tls_official_relaxed"
dirs <- list.dirs(BASE, recursive = FALSE)

for(d in dirs) {
  pdf_file <- file.path(d, "tls_score_official.pdf")
  png_file <- file.path(d, "tls_score_official.png")
  if(!file.exists(pdf_file)) next

  cat(basename(d), "... ")
  tryCatch({
    img <- image_read_pdf(pdf_file, density = 200)
    image_write(img, png_file, format = "png")
    cat("OK\n")
  }, error = function(e) cat("ERR:", conditionMessage(e), "\n"))
}
cat("Done\n")
