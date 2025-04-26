#' Bin Cells Along Pseudotime
#' @param seurat_obj Seurat object
#' @param pseudotime_col Metadata column name for pseudotime
#' @param n_bins Number of bins
#' @param bin_method Binning method ("quantile" or "equal_width")
#' @return Vector of bin assignments
#' @export
make_pt_bins <- function(seurat_obj, pseudotime_col, n_bins = 5, bin_method = "quantile") {
  pt <- seurat_obj@meta.data[[pseudotime_col]]
  bins <- rep(NA_integer_, length(pt))
  ok <- !is.na(pt)
  if (sum(ok) < 10) return(bins)
  if (bin_method == "equal_width") {
    breaks <- seq(min(pt[ok]), max(pt[ok]), length.out = n_bins + 1)
  } else {
    breaks <- quantile(pt[ok], probs = seq(0, 1, length.out = n_bins + 1))
  }
  bins[ok] <- cut(pt[ok], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bins
}

#' Gene-Pseudotime Correlation
#' @param seurat_obj Seurat object
#' @param pseudotime_col Metadata column name for pseudotime
#' @param assay Assay name
#' @param layer Layer name
#' @param method Correlation method
#' @return Data frame with gene correlations
#' @export
gene_pt_correlation <- function(seurat_obj, pseudotime_col, assay = "RNA", layer = "data", method = "spearman") {
  get_assay_data <- function(obj, assay, layer, slot) {
    if (packageVersion("SeuratObject") >= "5.0.0") {
      Seurat::GetAssayData(obj, assay = assay, layer = layer)
    } else {
      Seurat::GetAssayData(obj, assay = assay, slot = slot)
    }
  }
  expr <- get_assay_data(seurat_obj, assay, layer, layer)
  pt <- seurat_obj@meta.data[[pseudotime_col]]
  data.frame(Gene = rownames(expr),
             Cor = apply(expr, 1, cor, y = pt, method = method, use = "pairwise.complete.obs"),
             stringsAsFactors = FALSE)
}

#' Compute Edge-Pseudotime Correlations
#'
#' Computes correlations between edge weights and pseudotime for a given graph.
#'
#' @param graph An igraph object.
#' @param seurat_obj A Seurat object.
#' @param pseudotime_col Character, metadata column for pseudotime.
#' @param assay Character, assay to use.
#' @param layer Character, data layer to use.
#' @param method Character, correlation method.
#' @param edge_threshold Numeric, edge weight threshold.
#' @param n_cores Integer, number of cores for parallel computation.
#' @return A data frame with edge correlations.
#' @export
#' @importFrom stats cor sd
edge_pt_correlation <- function(graph, seurat_obj, pseudotime_col = "pseudotime",
                                assay = "RNA", layer = "data", method = "spearman",
                                edge_threshold = 0.95, n_cores = 1) {
  # Extract expression data and pseudotime
  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  pseudotime <- seurat_obj@meta.data[[pseudotime_col]]

  # Filter edges based on threshold
  edges <- igraph::as_data_frame(graph, what = "edges")
  edges <- edges[edges$weight >= quantile(edges$weight, edge_threshold), ]

  # Compute correlations
  cor_results <- parallel::mclapply(1:nrow(edges), function(i) {
    gene1 <- edges$from[i]
    gene2 <- edges$to[i]
    expr1 <- expr_data[gene1, ]
    expr2 <- expr_data[gene2, ]

    # Check for zero variance
    if (stats::sd(expr1) == 0 || stats::sd(expr2) == 0 || stats::sd(pseudotime) == 0) {
      return(data.frame(from = gene1, to = gene2, cor_pt = NA, stringsAsFactors = FALSE))
    }

    cor_val <- stats::cor(expr1 + expr2, pseudotime, method = method, use = "pairwise.complete.obs")
    data.frame(from = gene1, to = gene2, cor_pt = cor_val, stringsAsFactors = FALSE)
  }, mc.cores = n_cores)

  do.call(rbind, cor_results)
}

#' Edge Gain/Loss Analysis
#' @param bin_graphs List of graphs per bin
#' @return List of stable, gained, and lost edges
#' @export
edge_gain_loss <- function(bin_graphs) {
  edge_key <- function(g)
    if (is.null(g)) character(0) else
      apply(igraph::as_edgelist(g), 1, function(e) paste(sort(e), collapse = "_"))
  keys_by_bin <- lapply(bin_graphs, edge_key)
  list(stable = Reduce(intersect, keys_by_bin),
       gained = setdiff(keys_by_bin[[length(keys_by_bin)]], keys_by_bin[[1]]),
       lost = setdiff(keys_by_bin[[1]], keys_by_bin[[length(keys_by_bin)]]))
}

#' Compute Conditional MI Matrix in Parallel
#' @param expr_subset Expression data subset
#' @param pt_vec Pseudotime vector
#' @param n_cores Number of cores
#' @param nbins Number of bins
#' @return Conditional MI matrix
#' @export
compute_cmi_matrix_parallel <- function(expr_subset, pt_vec, n_cores = parallel::detectCores() - 10, nbins = 10) {
  disc_vec <- function(v, nbins = 10) {
    if (var(v, na.rm = TRUE) == 0) return(rep(0L, length(v)))
    nb <- min(nbins, length(unique(v)))
    infotheo::discretize(matrix(v, ncol = 1), disc = "equalfreq", nbins = nb)[, 1]
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
  rownames(cmi_mat) <- colnames(cmi_mat) <- rownames(expr_subset)
  cmi_mat
}

#' Compute Pseudotime-Aware MI Networks
#' @param seurat_obj Seurat object
#' @param cluster_id Cluster ID
#' @param cluster_field Metadata column name for clusters
#' @param pseudotime_col Metadata column name for pseudotime
#' @param assay_name Assay name (default: "RNA")
#' @param counts_layer Layer for counts data (default: "counts")
#' @param data_layer Layer for normalized data (default: "data")
#' @param n_bins Number of pseudotime bins
#' @param min_cells_bin Minimum cells per bin
#' @param top_n_genes Number of top variable genes
#' @param max_cells Maximum number of cells
#' @param min_expr_pct Minimum expression percentage
#' @param n_cores Number of cores
#' @param bin_method Binning method ("quantile" or "equal_width")
#' @return List of MI networks per pseudotime bin
#' @export
compute_pseudotime_mi_network <- function(seurat_obj, cluster_id, cluster_field, pseudotime_col,
                                          assay_name = "RNA", counts_layer = "counts", data_layer = "data",
                                          n_bins = 5, min_cells_bin = 75, top_n_genes = 5000,
                                          max_cells = 1500, min_expr_pct = 0.05, n_cores = parallel::detectCores() - 1,
                                          bin_method = "quantile") {
  if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
    stop("Pseudotime column '", pseudotime_col, "' not found in Seurat object metadata")
  }
  keep <- seurat_obj@meta.data[[cluster_field]] == cluster_id & !is.na(seurat_obj@meta.data[[pseudotime_col]])
  if (sum(keep) == 0) {
    stop("No cells found for cluster '", cluster_id, "' with non-NA pseudotime")
  }
  bins_full <- make_pt_bins(seurat_obj[, keep], pseudotime_col, n_bins, bin_method)

  bin_graphs <- lapply(1:n_bins, function(b) {
    cells <- colnames(seurat_obj)[keep & bins_full == b]
    if (length(cells) < min_cells_bin) {
      message("Bin ", b, ": skipped (", length(cells), " cells)")
      return(NULL)
    }
    expr <- get_expr_data(seurat_obj[, cells], cluster_id, cluster_field,
                          assay_name, counts_layer, data_layer,
                          top_n_genes, max_cells, min_expr_pct)
    mi <- compute_mi_matrix_parallel(expr, n_cores)
    gc()
    build_mi_graph(mi)
  })

  names(bin_graphs) <- paste0("PT", 1:n_bins)
  bin_graphs
}
