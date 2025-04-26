#' Compute Mutual Information Matrix in Parallel
#' @param expr_subset Expression data subset
#' @param n_cores Number of cores
#' @return MI matrix
#' @export
compute_mi_matrix_parallel <- function(expr_subset, n_cores = parallel::detectCores() - 10) {
  sd_vals <- apply(expr_subset, 1, function(x) sd(x, na.rm = TRUE))
  if (any(is.na(sd_vals) | sd_vals == 0)) {
    cat("Warning: Found genes with NA or zero standard deviation:\n")
    print(names(which(is.na(sd_vals) | sd_vals == 0)))
  }
  nbins_vec <- apply(expr_subset, 1, choose_adaptive_nbins)
  nbins <- round(median(nbins_vec, na.rm = TRUE))
  cat("Using adaptive nbins =", nbins, "\n")
  disc_expr <- infotheo::discretize(t(expr_subset), nbins = nbins)
  n_genes <- nrow(expr_subset)
  gene_pairs <- combn(n_genes, 2, simplify = FALSE)
  mi_values <- parallel::mclapply(gene_pairs, function(pair) {
    infotheo::mutinformation(disc_expr[, pair[1]], disc_expr[, pair[2]])
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  mi_mat <- matrix(0, n_genes, n_genes)
  idx <- combn(n_genes, 2)
  mi_mat[idx] <- mi_mat[t(idx)] <- unlist(mi_values)
  diag(mi_mat) <- 0
  rownames(mi_mat) <- colnames(mi_mat) <- rownames(expr_subset)
  mi_mat
}

#' Get Expression Data for a Cluster
#' @param seurat_obj Seurat object
#' @param cluster_id Cluster ID
#' @param cluster_field Metadata column name for clusters
#' @param assay_name Assay name (default: "RNA")
#' @param counts_layer Layer for counts data (default: "counts")
#' @param data_layer Layer for normalized data (default: "data")
#' @param top_n_genes Number of top variable genes
#' @param max_cells Maximum number of cells
#' @param min_expr_pct Minimum expression percentage
#' @return Expression matrix
#' @export
get_expr_data <- function(seurat_obj, cluster_id, cluster_field,
                          assay_name = "RNA", counts_layer = "counts", data_layer = "data",
                          top_n_genes = 5000, max_cells = 1500, min_expr_pct = 0.05) {
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop("Specified assay '", assay_name, "' not found in Seurat object")
  }
  if (!cluster_field %in% colnames(seurat_obj@meta.data)) {
    stop("Specified cluster field '", cluster_field, "' not found in Seurat object metadata")
  }
  cluster_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_field]] == cluster_id]
  if (length(cluster_cells) == 0) {
    stop("No cells found for cluster '", cluster_id, "' in field '", cluster_field, "'")
  }
  if (length(cluster_cells) > max_cells) {
    set.seed(123)
    cluster_cells <- sample(cluster_cells, max_cells)
  }

  # Use layer for Seurat v5, fall back to slot for v4
  get_assay_data <- function(obj, assay, layer, slot) {
    if (packageVersion("SeuratObject") >= "5.0.0") {
      Seurat::GetAssayData(obj, assay = assay, layer = layer)
    } else {
      Seurat::GetAssayData(obj, assay = assay, slot = slot)
    }
  }

  counts_data <- get_assay_data(seurat_obj[, cluster_cells], assay_name, counts_layer, counts_layer)
  if (!is.matrix(counts_data) && !inherits(counts_data, "Matrix")) {
    stop("Counts data is not a matrix. Check assay '", assay_name, "' and layer '", counts_layer, "'")
  }
  if (nrow(counts_data) < 2 || ncol(counts_data) < 2) {
    stop("Counts data has insufficient dimensions (", nrow(counts_data), " rows, ", ncol(counts_data), " columns)")
  }
  min_cells <- ceiling(min_expr_pct * ncol(counts_data))
  gene_expr_counts <- Matrix::rowSums(counts_data > 0)  # Explicitly use Matrix::rowSums for sparse matrices
  genes_to_keep <- names(gene_expr_counts[gene_expr_counts >= min_cells])
  if (length(genes_to_keep) < 2) {
    stop("Fewer than 2 genes remain after expression filter for cluster '", cluster_id, "'")
  }
  message("Cluster ", cluster_id, ": ", length(genes_to_keep), " genes kept after ", min_expr_pct * 100, "% filter")
  rm(counts_data); gc(verbose = FALSE)

  norm_data <- get_assay_data(seurat_obj[, cluster_cells], assay_name, data_layer, data_layer)
  if (!is.matrix(norm_data) && !inherits(norm_data, "Matrix")) {
    stop("Normalized data is not a matrix. Check assay '", assay_name, "' and layer '", data_layer, "'")
  }
  if (nrow(norm_data) < 2 || ncol(norm_data) < 2) {
    stop("Normalized data has insufficient dimensions (", nrow(norm_data), " rows, ", ncol(norm_data), " columns)")
  }
  norm_data <- norm_data[genes_to_keep, , drop = FALSE]
  gene_vars <- apply(norm_data, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n_genes, length(gene_vars))]
  filtered_genes <- filter_housekeeping(top_genes)
  if (length(filtered_genes) < 2) {
    stop("Fewer than 2 genes remain after housekeeping filter for cluster '", cluster_id, "'")
  }
  message("Cluster ", cluster_id, ": ", length(filtered_genes), " genes after housekeeping filter")
  as.matrix(norm_data[filtered_genes, , drop = FALSE])
}

#' Build MI Graph
#' @param mi_matrix Mutual information matrix
#' @param threshold_method Threshold method, either "percentile" or "fixed"
#' @param percentile Percentile for thresholding (used if threshold_method = "percentile")
#' @param fixed_threshold Fixed MI value for thresholding (used if threshold_method = "fixed")
#' @return igraph object
#' @export
build_mi_graph <- function(mi_matrix, threshold_method = "percentile",
                           percentile = 0.95, fixed_threshold = NULL) {
  if (!is.matrix(mi_matrix) || nrow(mi_matrix) != ncol(mi_matrix)) {
    stop("mi_matrix must be a square matrix")
  }
  if (!threshold_method %in% c("percentile", "fixed")) {
    stop("threshold_method must be 'percentile' or 'fixed'")
  }

  mi_vals <- mi_matrix[upper.tri(mi_matrix)]
  if (length(mi_vals) == 0) {
    stop("No valid MI values in the upper triangle of mi_matrix")
  }

  if (threshold_method == "percentile") {
    if (!is.numeric(percentile) || percentile <= 0 || percentile >= 1) {
      stop("percentile must be a number between 0 and 1")
    }
    threshold <- quantile(mi_vals, percentile, na.rm = TRUE)
  } else { # threshold_method == "fixed"
    if (is.null(fixed_threshold) || !is.numeric(fixed_threshold) || fixed_threshold < 0) {
      stop("fixed_threshold must be a non-negative numeric value when threshold_method = 'fixed'")
    }
    threshold <- fixed_threshold
  }

  adj_matrix <- matrix(0, nrow = nrow(mi_matrix), ncol = ncol(mi_matrix))
  adj_matrix[mi_matrix >= threshold] <- mi_matrix[mi_matrix >= threshold]
  diag(adj_matrix) <- 0
  rownames(adj_matrix) <- colnames(adj_matrix) <- rownames(mi_matrix)
  igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
}

#' Compute MI Network for a Cluster
#' @param seurat_obj Seurat object
#' @param cluster_id Cluster ID
#' @param cluster_field Metadata column name for clusters
#' @param assay_name Assay name (default: "RNA")
#' @param counts_layer Layer for counts data (default: "counts")
#' @param data_layer Layer for normalized data (default: "data")
#' @param top_n_genes Number of top variable genes
#' @param max_cells Maximum number of cells
#' @param min_expr_pct Minimum expression percentage
#' @param n_cores Number of cores
#' @return List containing graph, MI matrix, and genes
#' @export
compute_mi_network <- function(seurat_obj, cluster_id, cluster_field,
                               assay_name = "RNA", counts_layer = "counts", data_layer = "data",
                               top_n_genes = 5000, max_cells = 1500, min_expr_pct = 0.05,
                               n_cores = parallel::detectCores() - 1) {
  expr_subset <- get_expr_data(seurat_obj, cluster_id, cluster_field,
                               assay_name, counts_layer, data_layer,
                               top_n_genes, max_cells, min_expr_pct)
  mi_mat <- compute_mi_matrix_parallel(expr_subset, n_cores)
  g <- build_mi_graph(mi_mat, threshold_method = "percentile", percentile = 0.95)
  list(graph = g, mi_matrix = mi_mat, genes = rownames(expr_subset))
}

#' Percolate Graph
#'
#' Creates a consensus graph by filtering edges at the 95th percentile and retaining those in >=50% of iterations.
#'
#' @param graph Input igraph object.
#' @param n_iterations Number of percolation iterations.
#' @param edge_fraction Fraction of edges to retain per iteration.
#' @param freq_threshold Frequency threshold for edge retention.
#' @param percentile Percentile for initial edge filtering (default: 0.95).
#'
#' @return Consensus igraph object.
#'
#' @export
#' @importFrom igraph as_data_frame graph_from_data_frame E subgraph.edges list.edge.attributes
#' @importFrom dplyr filter group_by summarise
percolate_graph <- function(graph, n_iterations = 50, edge_fraction = 0.6, freq_threshold = 0.5, percentile = 0.95) {
  # Pre-filter edges at 95th percentile
  el_full <- igraph::as_data_frame(graph, what = "edges")
  top_thr <- quantile(el_full$weight, percentile, na.rm = TRUE)
  el_thresh <- el_full %>% dplyr::filter(weight >= top_thr)
  g_thresh <- igraph::graph_from_data_frame(el_thresh, directed = FALSE)

  # Percolation
  consensus_edge_list <- list()
  for (i in seq_len(n_iterations)) {
    set.seed(321 + i)
    edges_to_keep <- sample(igraph::E(g_thresh), size = round(length(igraph::E(g_thresh)) * edge_fraction))
    g_perc_iter <- igraph::subgraph.edges(g_thresh, edges_to_keep)
    consensus_edge_list[[i]] <- igraph::as_data_frame(g_perc_iter, what = "edges")
  }
  combined_edges <- dplyr::bind_rows(consensus_edge_list) %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(freq = dplyr::n(), mean_weight = mean(weight), .groups = "drop")
  consensus_edges <- combined_edges %>% dplyr::filter(freq >= freq_threshold * n_iterations)
  g_consensus <- igraph::graph_from_data_frame(consensus_edges, directed = FALSE)
  if ("mean_weight" %in% igraph::list.edge.attributes(g_consensus)) {
    igraph::E(g_consensus)$weight <- igraph::E(g_consensus)$mean_weight
  }
  g_consensus
}

#' Detect Communities
#' Detects communities using label propagation and neighborhood overlap.
#' @param graph Input igraph object.
#' @param min_size Minimum community size.
#' @param overlap_threshold Overlap threshold for adding genes (default: 0.1).
#' @return List with communities and their sizes.
#' @export
detect_communities <- function(graph, min_size = 10, overlap_threshold = 0.1) {
  set.seed(123)
  lc_comm <- cluster_label_prop(graph, weights = E(graph)$weight)
  lc_groups <- split(V(graph)$name, membership(lc_comm))
  lc_groups <- Filter(function(x) length(x) >= min_size, lc_groups)

  # Neighborhood overlap
  adj_matrix <- as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)
  for (node in V(graph)$name) {
    idx <- which(V(graph)$name == node)
    neighbors <- V(graph)$name[adj_matrix[idx, ] > 0]
    if (length(neighbors) == 0) next
    for (i in seq_along(lc_groups)) {
      comm_genes <- lc_groups[[i]]
      ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
      if (ov >= overlap_threshold && !(node %in% comm_genes)) {
        lc_groups[[i]] <- c(comm_genes, node)
      }
    }
  }
  large_communities <- Filter(function(x) length(x) >= min_size, lc_groups)
  list(communities = large_communities, sizes = sapply(large_communities, length))
}
