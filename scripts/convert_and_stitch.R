library(magick)

BASE <- "E:/GBM/results/tls_official_relaxed"
dirs <- list.dirs(BASE, recursive = FALSE)

# Step 1: Convert all PDFs to PNG
png_files <- c()
for(d in dirs) {
  pdf_file <- file.path(d, "tls_score_official.pdf")
  png_file <- file.path(d, "tls_score_official.png")
  if(!file.exists(pdf_file)) next

  cat(basename(d), "... ")
  tryCatch({
    img <- image_read_pdf(pdf_file, density = 200)
    image_write(img, png_file, format = "png")
    png_files <- c(png_files, png_file)
    cat("OK\n")
  }, error = function(e) cat("ERR:", conditionMessage(e), "\n"))
}

cat(sprintf("\nConverted %d PDFs to PNG\n", length(png_files)))

# Step 2: Stitch into 4x5 grids (20 per grid)
if(length(png_files) == 0) quit()

n_per_grid <- 20
n_grids <- ceiling(length(png_files) / n_per_grid)

for(g in seq_len(n_grids)) {
  start <- (g-1) * n_per_grid + 1
  end <- min(g * n_per_grid, length(png_files))
  batch <- png_files[start:end]

  cat(sprintf("Grid %d/%d: %d images...\n", g, n_grids, length(batch)))

  # Read and optionally resize all images in this batch
  imgs <- lapply(batch, function(f) {
    img <- image_read(f)
    # Resize to uniform size for grid
    image_scale(img, "600x600")
  })

  # Build 4-column grid
  # Fill rows: pad with blank if needed
  while(length(imgs) < n_per_grid) {
    imgs <- c(imgs, list(image_blank(600, 600, "white")))
  }

  # Arrange in rows of 4
  rows <- list()
  for(r in seq_len(5)) {
    row_start <- (r-1) * 4 + 1
    row_end <- r * 4
    rows[[r]] <- image_append(do.call(c, imgs[row_start:row_end]))
  }

  grid <- image_append(do.call(c, rows), stack = TRUE)

  out_file <- file.path(BASE, sprintf("tls_score_grid_%02d.png", g))
  image_write(grid, out_file, format = "png")
  cat(sprintf("  -> %s\n", out_file))
}

cat("\nDone!\n")
