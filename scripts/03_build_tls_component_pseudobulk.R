library(rhdf5)
library(Matrix)

root <- "E:/GBM"
h5_dir <- file.path(root, "spatial_data_visium", "spatial_data_visium", "anndata_consolidated")
c2l_dir <- file.path(root, "results", "c2l_core_v1")
tls_root <- file.path(root, "results", "tls_core")
out_root <- file.path(root, "results")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

build_tls_components <- function(
  tls_df,
  radius = 2,
  min_spots = 10,
  saddle_ratio = 0.85,
  edge_margin = 1
) {
  tls_only <- tls_df[tls_df[["TLS.region"]] == "TLS", , drop = FALSE]
  if (nrow(tls_only) < min_spots) return(NULL)

  xy <- cbind(
    as.numeric(tls_only[["array_col"]]),
    as.numeric(tls_only[["array_row"]]) * sqrt(3)
  )
  d <- as.matrix(stats::dist(xy))
  adj <- d <= radius + 1e-8
  diag(adj) <- TRUE
  tls_neighbors <- lapply(seq_len(nrow(tls_only)), function(i) which(adj[i, ]))
  tissue_df <- tls_df[tls_df[["in_tissue"]] == 1, , drop = FALSE]
  row_min <- min(as.numeric(tissue_df[["array_row"]]), na.rm = TRUE)
  row_max <- max(as.numeric(tissue_df[["array_row"]]), na.rm = TRUE)
  col_min <- min(as.numeric(tissue_df[["array_col"]]), na.rm = TRUE)
  col_max <- max(as.numeric(tissue_df[["array_col"]]), na.rm = TRUE)

  split_connected <- function(nodes, local_adj) {
    out <- list()
    seen <- rep(FALSE, length(nodes))
    map <- setNames(seq_along(nodes), nodes)
    for (ii in seq_along(nodes)) {
      if (seen[ii]) next
      queue <- ii
      seen[ii] <- TRUE
      keep <- integer(0)
      while (length(queue) > 0) {
        v <- queue[1]
        queue <- queue[-1]
        keep <- c(keep, nodes[v])
        nei_global <- intersect(local_adj[[nodes[v]]], nodes)
        nei <- unname(map[as.character(nei_global)])
        nei <- nei[!is.na(nei) & !seen[nei]]
        if (length(nei) > 0) {
          seen[nei] <- TRUE
          queue <- c(queue, nei)
        }
      }
      out[[length(out) + 1L]] <- sort(unique(keep))
    }
    out
  }

  absorb_small <- function(parts, local_adj, min_spots) {
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      sizes <- vapply(parts, length, integer(1))
      small <- which(sizes < min_spots)
      if (length(small) == 0 || length(parts) <= 1) break
      for (sid in small) {
        if (sid > length(parts)) next
        members <- parts[[sid]]
        border <- unique(unlist(lapply(members, function(m) setdiff(local_adj[[m]], members))))
        if (length(border) == 0) next
        hit <- which(vapply(parts, function(x) any(x %in% border), logical(1)))
        hit <- setdiff(hit, sid)
        if (length(hit) == 0) next
        target <- hit[[1]]
        parts[[target]] <- sort(unique(c(parts[[target]], members)))
        parts[[sid]] <- integer(0)
        changed <- TRUE
      }
      parts <- Filter(length, parts)
    }
    parts
  }

  find_root <- function(parent, x) {
    while (parent[x] != x) x <- parent[x]
    x
  }

  union_root <- function(parent, a, b) {
    ra <- find_root(parent, a)
    rb <- find_root(parent, b)
    if (ra != rb) parent[rb] <- ra
    parent
  }

  is_bad_edge_line <- function(nodes, local_adj, local_df, local_xy) {
    if (length(nodes) < min_spots) return(TRUE)
    sub <- local_df[nodes, , drop = FALSE]
    rr <- as.numeric(sub[["array_row"]])
    cc <- as.numeric(sub[["array_col"]])
    touch_edge <- any(
      rr <= (row_min + edge_margin) | rr >= (row_max - edge_margin) |
        cc <= (col_min + edge_margin) | cc >= (col_max - edge_margin)
    )
    neigh_deg <- vapply(nodes, function(m) {
      sum(local_adj[[m]] %in% nodes) - 1L
    }, integer(1))
    row_span <- max(rr) - min(rr) + 1
    col_span <- max(cc) - min(cc) + 1
    occupancy <- length(nodes) / max(1, row_span * col_span)
    eig <- eigen(stats::cov(local_xy[nodes, , drop = FALSE]), symmetric = TRUE, only.values = TRUE)$values
    elong <- if (sum(eig, na.rm = TRUE) > 0) max(eig, na.rm = TRUE) / sum(eig, na.rm = TRUE) else 1
    touch_edge && max(neigh_deg) <= 2 && (occupancy < 0.45 || elong > 0.93)
  }

  feature_cols <- intersect(
    c("TLS.score", "Plasma_B.cells_T.cells", "LC.50sig", "Plasma_B_cells", "T_cells"),
    colnames(tls_only)
  )
  if (length(feature_cols) == 0) feature_cols <- "TLS.score"
  feat <- as.data.frame(lapply(feature_cols, function(col) {
    x <- suppressWarnings(as.numeric(tls_only[[col]]))
    x[!is.finite(x)] <- 0
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }))
  colnames(feat) <- feature_cols
  raw_strength <- rowMeans(feat, na.rm = TRUE)
  strength <- vapply(seq_len(nrow(tls_only)), function(i) {
    mean(raw_strength[tls_neighbors[[i]]], na.rm = TRUE)
  }, numeric(1))

  connected <- split_connected(seq_len(nrow(tls_only)), tls_neighbors)
  parts <- list()

  for (cc_nodes in connected) {
    if (length(cc_nodes) < min_spots) next
    local_map <- setNames(seq_along(cc_nodes), cc_nodes)
    local_strength <- strength[cc_nodes]
    local_raw <- raw_strength[cc_nodes]

    peak_mask <- vapply(cc_nodes, function(node) {
      nei <- setdiff(intersect(tls_neighbors[[node]], cc_nodes), node)
      length(nei) == 0 || all(strength[node] >= strength[nei] - 1e-8)
    }, logical(1))
    peak_nodes <- cc_nodes[peak_mask]
    if (length(peak_nodes) == 0) peak_nodes <- cc_nodes[which.max(local_strength)]

    peak_parent <- seq_along(peak_nodes)
    names(peak_parent) <- peak_nodes
    peak_adj_idx <- which(adj[peak_nodes, peak_nodes, drop = FALSE], arr.ind = TRUE)
    peak_adj_idx <- peak_adj_idx[peak_adj_idx[, 1] < peak_adj_idx[, 2], , drop = FALSE]
    if (nrow(peak_adj_idx) > 0) {
      for (k in seq_len(nrow(peak_adj_idx))) {
        a <- peak_adj_idx[k, 1]
        b <- peak_adj_idx[k, 2]
        if (abs(strength[peak_nodes[a]] - strength[peak_nodes[b]]) <= 0.02) {
          peak_parent <- union_root(peak_parent, a, b)
        }
      }
    }
    peak_group <- vapply(seq_along(peak_nodes), function(i) find_root(peak_parent, i), integer(1))
    peak_group <- match(peak_group, sort(unique(peak_group)))

    basin_label <- rep(NA_integer_, length(cc_nodes))
    names(basin_label) <- cc_nodes
    peak_lookup <- rep(NA_integer_, nrow(tls_only))
    peak_lookup[peak_nodes] <- peak_group

    for (node in cc_nodes[order(local_strength, local_raw, decreasing = TRUE)]) {
      if (!is.na(peak_lookup[node])) {
        basin_label[as.character(node)] <- peak_lookup[node]
        next
      }
      cur <- node
      seen <- integer(0)
      repeat {
        seen <- c(seen, cur)
        nei <- intersect(tls_neighbors[[cur]], cc_nodes)
        better <- nei[strength[nei] > strength[cur] + 1e-8]
        if (length(better) == 0) break
        ord <- order(strength[better], raw_strength[better], decreasing = TRUE)
        nxt <- better[ord[[1]]]
        if (nxt %in% seen) break
        cur <- nxt
      }
      if (!is.na(peak_lookup[cur])) {
        basin_label[as.character(node)] <- peak_lookup[cur]
      } else {
        best_peak <- which.max(strength[peak_nodes] - d[node, peak_nodes] * 1e-6)
        basin_label[as.character(node)] <- peak_group[best_peak]
      }
    }

    basin_parts <- split(cc_nodes, basin_label[as.character(cc_nodes)])
    basin_parts <- lapply(basin_parts, sort)
    basin_parts <- Filter(length, basin_parts)
    if (length(basin_parts) > 1) {
      basin_parent <- seq_along(basin_parts)
      edge_idx <- which(adj[cc_nodes, cc_nodes, drop = FALSE], arr.ind = TRUE)
      edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]
      if (nrow(edge_idx) > 0) {
        edge_df <- data.frame(
          u = cc_nodes[edge_idx[, 1]],
          v = cc_nodes[edge_idx[, 2]]
        )
        edge_df$a <- match(basin_label[as.character(edge_df$u)], names(basin_parts))
        edge_df$b <- match(basin_label[as.character(edge_df$v)], names(basin_parts))
        edge_df <- edge_df[edge_df$a != edge_df$b, , drop = FALSE]
        if (nrow(edge_df) > 0) {
          peak_strength <- vapply(basin_parts, function(nodes) max(strength[nodes], na.rm = TRUE), numeric(1))
          pair_key <- paste(pmin(edge_df$a, edge_df$b), pmax(edge_df$a, edge_df$b), sep = "__")
          for (pk in unique(pair_key)) {
            sub <- edge_df[pair_key == pk, , drop = FALSE]
            a <- min(sub$a)
            b <- max(sub$b)
            saddle <- max(pmin(strength[sub$u], strength[sub$v]), na.rm = TRUE)
            if (is.finite(saddle) && saddle >= saddle_ratio * min(peak_strength[a], peak_strength[b])) {
              basin_parent <- union_root(basin_parent, a, b)
            }
          }
        }
      }
      basin_group <- vapply(seq_along(basin_parts), function(i) find_root(basin_parent, i), integer(1))
      basin_group <- match(basin_group, sort(unique(basin_group)))
      basin_parts <- split(unlist(basin_parts, use.names = FALSE), basin_group[rep(seq_along(basin_parts), lengths(basin_parts))])
      basin_parts <- lapply(basin_parts, function(x) sort(unique(x)))
    }

    basin_parts <- absorb_small(basin_parts, tls_neighbors, min_spots)
    basin_parts <- Filter(function(z) !is_bad_edge_line(z, tls_neighbors, tls_only, xy), basin_parts)
    basin_parts <- Filter(function(z) length(z) >= min_spots, basin_parts)
    parts <- c(parts, basin_parts)
  }

  if (length(parts) == 0) return(NULL)

  out <- do.call(rbind, lapply(seq_along(parts), function(i) {
    sub <- tls_only[parts[[i]], , drop = FALSE]
    sub$component_id <- i - 1L
    sub
  }))
  out
}

read_h5_counts <- function(path) {
  genes <- h5read(path, "var/_index")
  barcodes <- h5read(path, "obs/_index")
  xdata <- as.numeric(h5read(path, "X/data"))
  idx <- as.integer(h5read(path, "X/indices")) + 1L
  indptr <- as.integer(h5read(path, "X/indptr"))
  row_id <- rep(seq_len(length(barcodes)), diff(indptr))

  mat_obs_gene <- sparseMatrix(
    i = row_id,
    j = idx,
    x = xdata,
    dims = c(length(barcodes), length(genes))
  )
  rownames(mat_obs_gene) <- barcodes
  colnames(mat_obs_gene) <- genes

  list(
    counts = t(mat_obs_gene),
    genes = genes,
    barcodes = barcodes
  )
}

read_c2l_mean <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}

sample_ids <- intersect(
  intersect(
    sub("\\.h5ad$", "", basename(list.files(h5_dir, pattern = "\\.h5ad$", full.names = FALSE))),
    basename(list.dirs(tls_root, recursive = FALSE, full.names = FALSE))
  ),
  basename(list.dirs(c2l_dir, recursive = FALSE, full.names = FALSE))
)
sample_ids <- sort(sample_ids)

component_meta <- list()
component_spots <- list()
count_list <- list()
comp_list <- list()
all_ct <- NULL
all_genes <- character(0)

for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[[idx]]
  tls_csv <- file.path(tls_root, sid, "tls_spot_scores_official_relaxed.csv")
  h5ad <- file.path(h5_dir, paste0(sid, ".h5ad"))
  c2l_csv <- file.path(c2l_dir, sid, "cell2loc_mean.csv")
  if (!file.exists(tls_csv) || !file.exists(h5ad) || !file.exists(c2l_csv)) next

  tls_df <- read.csv(tls_csv, check.names = FALSE)
  comp_df <- build_tls_components(tls_df, radius = 2, min_spots = 10)
  if (is.null(comp_df)) next

  dat <- read_h5_counts(h5ad)
  c2l_mat <- read_c2l_mean(c2l_csv)
  if (is.null(all_ct)) all_ct <- colnames(c2l_mat)
  all_genes <- union(all_genes, rownames(dat$counts))

  shared <- Reduce(intersect, list(comp_df[["barcode"]], colnames(dat$counts), rownames(c2l_mat)))
  comp_df <- comp_df[comp_df[["barcode"]] %in% shared, , drop = FALSE]
  if (nrow(comp_df) == 0) next

  split_comp <- split(comp_df, comp_df[["component_id"]])
  kept <- 0L

  for (cid_name in names(split_comp)) {
    sub <- split_comp[[cid_name]]
    cid <- as.integer(cid_name)
    bc <- sub[["barcode"]]
    bc <- bc[bc %in% colnames(dat$counts)]
    if (length(bc) < 10) next

    unit_id <- paste0(sid, "__c", cid)
    count_vec <- Matrix::rowSums(dat$counts[, bc, drop = FALSE])
    count_vec <- setNames(as.numeric(count_vec), rownames(dat$counts))
    count_list[[unit_id]] <- count_vec
    comp_vec <- colMeans(c2l_mat[bc, , drop = FALSE])
    comp_vec <- comp_vec[all_ct]
    comp_list[[unit_id]] <- comp_vec
    kept <- kept + 1L

    component_meta[[unit_id]] <- data.frame(
      unit_id = unit_id,
      sample = sid,
      component_id = cid,
      n_spots = length(bc),
      tls_score_mean = mean(sub[["TLS.score"]], na.rm = TRUE),
      array_col_mean = mean(as.numeric(sub[["array_col"]]), na.rm = TRUE),
      array_row_mean = mean(as.numeric(sub[["array_row"]]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    component_spots[[unit_id]] <- data.frame(
      unit_id = unit_id,
      sample = sid,
      component_id = cid,
      barcode = bc,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("[%03d/%03d] %s: %d components kept\n", idx, length(sample_ids), sid, kept))
}

if (length(component_meta) == 0 || length(count_list) == 0 || length(comp_list) == 0) {
  stop(sprintf("No TLS components built. Checked tls_root=%s and h5_dir=%s", tls_root, h5_dir))
}

meta_df <- do.call(rbind, component_meta)
spots_df <- do.call(rbind, component_spots)
count_mat <- do.call(cbind, lapply(count_list, function(x) {
  out <- setNames(numeric(length(all_genes)), all_genes)
  out[names(x)] <- x
  out
}))
comp_mat <- do.call(rbind, lapply(comp_list, function(x) x[all_ct]))
colnames(comp_mat) <- all_ct

saveRDS(count_mat, file.path(out_root, "tls_pseudobulk_counts_by_component.rds"))
saveRDS(comp_mat, file.path(out_root, "tls_pseudobulk_c2l_by_component.rds"))
write.csv(meta_df, file.path(out_root, "tls_pseudobulk_component_metadata.csv"), row.names = FALSE)
write.csv(spots_df, file.path(out_root, "tls_component_spot_map.csv"), row.names = FALSE)

summary_df <- data.frame(
  n_units = nrow(meta_df),
  n_samples = length(unique(meta_df$sample)),
  median_spots = median(meta_df$n_spots),
  mean_spots = mean(meta_df$n_spots),
  min_spots = min(meta_df$n_spots),
  max_spots = max(meta_df$n_spots)
)
write.csv(summary_df, file.path(out_root, "tls_pseudobulk_component_summary.csv"), row.names = FALSE)

cat("Saved TLS pseudobulk component outputs\n")
print(summary_df)
