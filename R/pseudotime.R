#' Bin Cells Along Pseudotime
#' @param seurat_obj Seurat object
#' @param pseudotime_col Metadata column name for pseudotime
#' @param n_bins Number of bins
#' @param bin_method Binning method ("quantile" or "equal_width")
#' @return Integer vector of bin assignments (length = ncol(seurat_obj)); NA for cells without pseudotime
#' @export
make_pt_bins <- function(seurat_obj, pseudotime_col, n_bins = 5, bin_method = "quantile") {
  pt <- seurat_obj@meta.data[[pseudotime_col]]
  bins <- rep(NA_integer_, length(pt))
  ok <- !is.na(pt)

  if (sum(ok) < 10) return(bins)
  if (!is.numeric(pt[ok])) pt[ok] <- as.numeric(pt[ok])

  if (bin_method == "equal_width") {
    breaks <- seq(min(pt[ok]), max(pt[ok]), length.out = n_bins + 1)
  } else {
    breaks <- stats::quantile(pt[ok], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE, type = 7)
  }

  # guard against duplicated breaks (ties in pseudotime)
  breaks <- unique(as.numeric(breaks))
  if (length(breaks) < 3) {
    # not enough unique values to bin
    return(bins)
  }

  bins[ok] <- as.integer(cut(pt[ok], breaks = breaks, include.lowest = TRUE, labels = FALSE))
  bins
}

# internal helper: Seurat v4 slot vs v5 layer
.get_assay_data <- function(obj, assay, layer_or_slot) {
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      utils::packageVersion("SeuratObject") >= "5.0.0") {
    Seurat::GetAssayData(obj, assay = assay, layer = layer_or_slot)
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = layer_or_slot)
  }
}

#' Gene-Pseudotime Correlation
#' @param seurat_obj Seurat object
#' @param pseudotime_col Metadata column name for pseudotime
#' @param assay Assay name
#' @param layer Layer (Seurat v5) or slot (Seurat v4)
#' @param method Correlation method
#' @return Data frame with gene correlations
#' @export
#' @importFrom stats cor
gene_pt_correlation <- function(seurat_obj, pseudotime_col,
                                assay = "RNA", layer = "data", method = "spearman") {
  expr <- .get_assay_data(seurat_obj, assay, layer)
  pt <- seurat_obj@meta.data[[pseudotime_col]]

  data.frame(
    Gene = rownames(expr),
    Cor  = apply(expr, 1, stats::cor, y = pt, method = method, use = "pairwise.complete.obs"),
    stringsAsFactors = FALSE
  )
}

#' Compute Edge-Pseudotime Correlations
#'
#' Correlates (expr(gene1) + expr(gene2)) with pseudotime for top-weight edges.
#'
#' @param graph igraph object (must have edge attribute 'weight' or will use 1)
#' @param seurat_obj Seurat object
#' @param pseudotime_col Metadata column for pseudotime
#' @param assay Assay name
#' @param layer Layer (v5) / slot (v4)
#' @param method Correlation method
#' @param edge_threshold Quantile threshold for edge weights (0-1)
#' @param n_cores Cores for mclapply
#' @return Data frame with edge correlations
#' @export
#' @importFrom stats cor sd quantile
#' @importFrom igraph as_data_frame edge_attr_names E
edge_pt_correlation <- function(graph, seurat_obj, pseudotime_col = "pseudotime",
                                assay = "RNA", layer = "data", method = "spearman",
                                edge_threshold = 0.95, n_cores = 1) {
  expr_data <- .get_assay_data(seurat_obj, assay, layer)
  pseudotime <- seurat_obj@meta.data[[pseudotime_col]]

  edges <- igraph::as_data_frame(graph, what = "edges")
  if (!("weight" %in% colnames(edges))) {
    edges$weight <- 1
  }
  edges <- edges[is.finite(edges$weight) & !is.na(edges$weight), , drop = FALSE]
  if (nrow(edges) == 0) {
    return(data.frame(from = character(), to = character(), cor_pt = numeric(), stringsAsFactors = FALSE))
  }

  thr <- stats::quantile(edges$weight, probs = edge_threshold, na.rm = TRUE)
  edges <- edges[edges$weight >= thr, , drop = FALSE]
  if (nrow(edges) == 0) {
    return(data.frame(from = character(), to = character(), cor_pt = numeric(), stringsAsFactors = FALSE))
  }

  # only edges whose genes exist in expression matrix
  keep <- edges$from %in% rownames(expr_data) & edges$to %in% rownames(expr_data)
  edges <- edges[keep, , drop = FALSE]
  if (nrow(edges) == 0) {
    return(data.frame(from = character(), to = character(), cor_pt = numeric(), stringsAsFactors = FALSE))
  }

  # optional: restrict to cells with pseudotime
  ok_cells <- !is.na(pseudotime)
  if (sum(ok_cells) < 10) {
    return(data.frame(from = edges$from, to = edges$to, cor_pt = NA_real_, stringsAsFactors = FALSE))
  }
  expr_data <- expr_data[, ok_cells, drop = FALSE]
  pseudotime <- pseudotime[ok_cells]

  cor_results <- parallel::mclapply(seq_len(nrow(edges)), function(i) {
    gene1 <- edges$from[i]
    gene2 <- edges$to[i]
    expr1 <- as.numeric(expr_data[gene1, ])
    expr2 <- as.numeric(expr_data[gene2, ])

    if (stats::sd(expr1, na.rm = TRUE) == 0 ||
        stats::sd(expr2, na.rm = TRUE) == 0 ||
        stats::sd(pseudotime, na.rm = TRUE) == 0) {
      return(data.frame(from = gene1, to = gene2, cor_pt = NA_real_, stringsAsFactors = FALSE))
    }

    cor_val <- stats::cor(expr1 + expr2, pseudotime, method = method, use = "pairwise.complete.obs")
    data.frame(from = gene1, to = gene2, cor_pt = cor_val, stringsAsFactors = FALSE)
  }, mc.cores = n_cores)

  do.call(rbind, cor_results)
}

#' Edge Gain/Loss Analysis
#' @param bin_graphs List of graphs per bin (NULL allowed)
#' @return List of stable, gained, lost edge keys
#' @export
edge_gain_loss <- function(bin_graphs) {
  edge_key <- function(g) {
    if (is.null(g)) return(character(0))
    el <- igraph::as_edgelist(g, names = TRUE)
    apply(el, 1, function(e) paste(sort(e), collapse = "_"))
  }
  keys_by_bin <- lapply(bin_graphs, edge_key)
  list(
    stable = Reduce(intersect, keys_by_bin),
    gained = setdiff(keys_by_bin[[length(keys_by_bin)]], keys_by_bin[[1]]),
    lost   = setdiff(keys_by_bin[[1]], keys_by_bin[[length(keys_by_bin)]])
  )
}

#' Compute Conditional MI Matrix in Parallel
#' @param expr_subset Expression matrix (genes x cells)
#' @param pt_vec Pseudotime vector (length = cells)
#' @param n_cores Cores
#' @param nbins Discretization bins
#' @return Conditional MI matrix
#' @export
#' @importFrom stats var
compute_cmi_matrix_parallel <- function(expr_subset, pt_vec,
                                        n_cores = parallel::detectCores() - 10,
                                        nbins = 10) {
  if (!requireNamespace("infotheo", quietly = TRUE)) stop("Package 'infotheo' is required.")
  if (!requireNamespace("entropy", quietly = TRUE)) stop("Package 'entropy' is required.")

  disc_vec <- function(v, nbins = 10) {
    vv <- as.numeric(v)
    if (all(is.na(vv))) return(rep(0L, length(vv)))
    vv[is.na(vv)] <- stats::median(vv, na.rm = TRUE)
    if (stats::var(vv, na.rm = TRUE) == 0) return(rep(0L, length(vv)))
    nb <- min(nbins, length(unique(vv)))
    infotheo::discretize(matrix(vv, ncol = 1), disc = "equalfreq", nbins = nb)[, 1]
  }

  z_disc <- disc_vec(pt_vec, nbins)
  n_genes <- nrow(expr_subset)
  pairs <- combn(n_genes, 2, simplify = FALSE)

  cmi_vals <- parallel::mclapply(pairs, function(p) {
    xi <- disc_vec(expr_subset[p[1], ], nbins)
    yi <- disc_vec(expr_subset[p[2], ], nbins)
    entropy::mi.empirical(xi, yi, z_disc)
  }, mc.cores = n_cores, mc.set.seed = TRUE)

  cmi_mat <- matrix(0, n_genes, n_genes)
  idx <- combn(n_genes, 2)
  cmi_mat[idx] <- cmi_mat[t(idx)] <- unlist(cmi_vals)
  diag(cmi_mat) <- 0
  rownames(cmi_mat) <- colnames(cmi_mat) <- rownames(expr_subset)
  cmi_mat
}

#' Compute Pseudotime-Aware MI Networks (per pseudotime bin)
#' @param seurat_obj Seurat object
#' @param cluster_id Cluster ID value
#' @param cluster_field Metadata column name for clusters
#' @param pseudotime_col Metadata column name for pseudotime
#' @param assay_name Assay name (default: "RNA")
#' @param counts_layer Layer/slot for counts (default: "counts")
#' @param data_layer Layer/slot for normalized data (default: "data")
#' @param n_bins Number of pseudotime bins
#' @param min_cells_bin Minimum cells per bin
#' @param top_n_genes Top variable genes
#' @param max_cells Max cells (subsampling inside get_expr_data)
#' @param min_expr_pct Minimum expression percentage
#' @param n_cores Cores
#' @param bin_method "quantile" or "equal_width"
#' @return Named list of igraph objects per bin (NULL where skipped)
#' @export
compute_pseudotime_mi_network <- function(seurat_obj, cluster_id, cluster_field, pseudotime_col,
                                         assay_name = "RNA", counts_layer = "counts", data_layer = "data",
                                         n_bins = 5, min_cells_bin = 75, top_n_genes = 5000,
                                         max_cells = 1500, min_expr_pct = 0.05,
                                         n_cores = parallel::detectCores() - 1,
                                         bin_method = "quantile") {
  if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
    stop("Pseudotime column '", pseudotime_col, "' not found in Seurat object metadata")
  }
  if (!cluster_field %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster field '", cluster_field, "' not found in Seurat object metadata")
  }

  keep <- seurat_obj@meta.data[[cluster_field]] == cluster_id &
    !is.na(seurat_obj@meta.data[[pseudotime_col]])

  if (sum(keep) == 0) {
    stop("No cells found for cluster '", cluster_id, "' with non-NA pseudotime")
  }

  # bins aligned to full object length (fixes recycling bug)
  bins_all <- rep(NA_integer_, ncol(seurat_obj))
  bins_all[keep] <- make_pt_bins(seurat_obj[, keep, drop = FALSE], pseudotime_col, n_bins, bin_method)

  bin_graphs <- lapply(seq_len(n_bins), function(b) {
    cells <- colnames(seurat_obj)[which(bins_all == b)]
    if (length(cells) < min_cells_bin) {
      message("Bin ", b, ": skipped (", length(cells), " cells)")
      return(NULL)
    }
    expr <- get_expr_data(
      seurat_obj[, cells, drop = FALSE],
      cluster_id = cluster_id,
      cluster_field = cluster_field,
      assay_name = assay_name,
      counts_layer = counts_layer,
      data_layer = data_layer,
      top_n_genes = top_n_genes,
      max_cells = max_cells,
      min_expr_pct = min_expr_pct
    )
    mi <- compute_mi_matrix_parallel(expr, n_cores = n_cores)
    gc()
    build_mi_graph(mi, threshold_method = "percentile", percentile = 0.95)
  })

  names(bin_graphs) <- paste0("PT", seq_len(n_bins))
  bin_graphs
}

