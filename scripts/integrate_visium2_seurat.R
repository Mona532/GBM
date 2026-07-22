library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_visium2_seurat.R <input_dir> <output_dir> [sample_regex]")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
sample_regex <- if (length(args) >= 3) args[[3]] else "^#"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

region_map <- c(
  C = "Cortex",
  T = "Tumor",
  TC = "TumorCore",
  TI = "TumorInfiltration"
)

parse_sample_metadata <- function(sample_name) {
  clean <- sub("^#", "", sample_name)
  parts <- strsplit(clean, "_", fixed = TRUE)[[1]]
  patient <- parts[1]
  region_code <- parts[length(parts) - 1]
  idh_status <- if (grepl("IDHMutant", clean, ignore.case = TRUE)) "IDHmutant" else "IDHunknown"
  region_label <- if (region_code %in% names(region_map)) region_map[[region_code]] else region_code
  data.frame(
    sample_id = clean,
    sample_dir_name = sample_name,
    patient_id = patient,
    region_code = region_code,
    region_label = region_label,
    idh_status = idh_status,
    cohort = "UniversityClinicFreiburg",
    stringsAsFactors = FALSE
  )
}

load_one_sample <- function(sample_dir) {
  sample_name <- basename(sample_dir)
  outs_dir <- file.path(sample_dir, "outs")
  h5_path <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) {
    stop("Missing filtered_feature_bc_matrix.h5 for ", sample_name)
  }

  meta <- parse_sample_metadata(sample_name)
  slice_name <- meta$sample_id

  obj <- Load10X_Spatial(
    data.dir = outs_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = slice_name,
    filter.matrix = TRUE
  )

  obj <- RenameCells(obj, add.cell.id = slice_name)
  obj$sample_id <- meta$sample_id
  obj$sample_dir_name <- meta$sample_dir_name
  obj$patient_id <- meta$patient_id
  obj$region_code <- meta$region_code
  obj$region_label <- meta$region_label
  obj$idh_status <- meta$idh_status
  obj$cohort <- meta$cohort

  list(
    object = obj,
    summary = data.frame(
      sample_id = meta$sample_id,
      sample_dir_name = meta$sample_dir_name,
      patient_id = meta$patient_id,
      region_code = meta$region_code,
      region_label = meta$region_label,
      idh_status = meta$idh_status,
      n_spots = ncol(obj),
      n_features = nrow(obj),
      image_names = paste(Images(obj), collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

sample_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl(sample_regex, basename(sample_dirs))]
sample_dirs <- sort(sample_dirs)

if (!length(sample_dirs)) {
  stop("No sample directories matched regex: ", sample_regex)
}

loaded <- lapply(sample_dirs, load_one_sample)
objects <- lapply(loaded, `[[`, "object")
summary_df <- do.call(rbind, lapply(loaded, `[[`, "summary"))

merged <- objects[[1]]
if (length(objects) > 1) {
  merged <- merge(
    x = objects[[1]],
    y = objects[2:length(objects)],
    merge.data = FALSE
  )
}

merged@misc$integration_info <- list(
  source_dir = input_dir,
  n_samples = length(objects),
  sample_ids = summary_df$sample_id
)

saveRDS(summary_df, file.path(output_dir, "visium2_sample_summary.rds"))
write.csv(summary_df, file.path(output_dir, "visium2_sample_summary.csv"), row.names = FALSE)
saveRDS(merged, file.path(output_dir, "visium2_merged_raw_seurat.rds"))

message("[done] samples=", nrow(summary_df), " spots=", ncol(merged), " genes=", nrow(merged))
