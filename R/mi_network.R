#' Compute Mutual Information Matrix in Parallel
#' @param expr_subset Expression data subset (genes x cells)
#' @param n_cores Number of cores
#' @param chunk_size Number of gene pairs per chunk (default 200000)
#' @return MI matrix
#' @export
compute_mi_matrix_parallel <- function(expr_subset,
                                       n_cores = max(1L, parallel::detectCores() - 1L),
                                       chunk_size = 200000L) {

  if (!is.matrix(expr_subset)) expr_subset <- as.matrix(expr_subset)
  if (nrow(expr_subset) < 2 || ncol(expr_subset) < 2) {
    stop("expr_subset must be genes x cells with at least 2 genes and 2 cells.")
  }

  nbins_vec <- apply(expr_subset, 1, choose_adaptive_nbins)
  nbins <- as.integer(round(stats::median(nbins_vec, na.rm = TRUE)))
  if (!is.finite(nbins) || nbins < 2) nbins <- 5L

  disc_expr <- infotheo::discretize(t(expr_subset), nbins = nbins)

  n_genes <- nrow(expr_subset)
  idx <- utils::combn(n_genes, 2)
  n_pairs <- ncol(idx)

  chunks <- split(seq_len(n_pairs), ceiling(seq_len(n_pairs) / chunk_size))

  # PSOCK cluster is portable (works on HPC + Windows)
  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterSetRNGStream(cl, 123)
  parallel::clusterExport(cl, c("disc_expr", "idx"), envir = environment())
  parallel::clusterEvalQ(cl, library(infotheo))

  mi_vals <- unlist(parallel::parLapply(cl, chunks, function(cols) {
    vapply(cols, function(j) {
      infotheo::mutinformation(disc_expr[, idx[1, j]],
                               disc_expr[, idx[2, j]])
    }, numeric(1))
  }), use.names = FALSE)

  mi_mat <- matrix(0, n_genes, n_genes)
  mi_mat[idx] <- mi_vals
  mi_mat[t(idx)] <- mi_vals
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
#' @return Expression matrix (genes x cells)
#' @export
get_expr_data <- function(seurat_obj, cluster_id, cluster_field,
                          assay_name = "RNA",
                          counts_layer = "counts",
                          data_layer = "data",
                          top_n_genes = 5000,
                          max_cells = 1500,
                          min_expr_pct = 0.05) {

  if (!assay_name %in% names(seurat_obj@assays)) {
    stop("Specified assay '", assay_name, "' not found in Seurat object.")
  }
  if (!cluster_field %in% colnames(seurat_obj@meta.data)) {
    stop("Specified cluster field '", cluster_field, "' not found in Seurat metadata.")
  }

  cluster_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_field]] == cluster_id]
  if (length(cluster_cells) == 0) {
    stop("No cells found for cluster '", cluster_id, "' in field '", cluster_field, "'.")
  }
  if (length(cluster_cells) > max_cells) {
    set.seed(123)
    cluster_cells <- sample(cluster_cells, max_cells)
  }

  get_assay_data <- function(obj, assay, layer, slot) {
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      Seurat::GetAssayData(obj, assay = assay, layer = layer)
    } else {
      Seurat::GetAssayData(obj, assay = assay, slot = slot)
    }
  }

  counts_data <- get_assay_data(seurat_obj[, cluster_cells], assay_name, counts_layer, counts_layer)

  min_cells <- ceiling(min_expr_pct * ncol(counts_data))
  gene_expr_counts <- Matrix::rowSums(counts_data > 0)
  genes_to_keep <- names(gene_expr_counts[gene_expr_counts >= min_cells])
  if (length(genes_to_keep) < 2) {
    stop("Fewer than 2 genes remain after expression filter for cluster '", cluster_id, "'.")
  }

  norm_data <- get_assay_data(seurat_obj[, cluster_cells], assay_name, data_layer, data_layer)
  norm_data <- norm_data[genes_to_keep, , drop = FALSE]

  gene_vars <- apply(norm_data, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n_genes, length(gene_vars))]
  filtered_genes <- filter_housekeeping(top_genes)
  if (length(filtered_genes) < 2) {
    stop("Fewer than 2 genes remain after housekeeping filter for cluster '", cluster_id, "'.")
  }

  as.matrix(norm_data[filtered_genes, , drop = FALSE])
}

#' Build MI Graph
#' @param mi_matrix Mutual information matrix
#' @param threshold_method Threshold method, either "percentile" or "fixed"
#' @param percentile Percentile for thresholding (used if threshold_method = "percentile")
#' @param fixed_threshold Fixed MI value for thresholding (used if threshold_method = "fixed")
#' @return igraph object
#' @export
build_mi_graph <- function(mi_matrix,
                           threshold_method = "percentile",
                           percentile = 0.95,
                           fixed_threshold = NULL) {

  if (!is.matrix(mi_matrix) || nrow(mi_matrix) != ncol(mi_matrix)) {
    stop("mi_matrix must be a square matrix.")
  }
  if (!threshold_method %in% c("percentile", "fixed")) {
    stop("threshold_method must be 'percentile' or 'fixed'.")
  }

  mi_vals <- mi_matrix[upper.tri(mi_matrix)]
  if (length(mi_vals) == 0) stop("No MI values found (upper triangle empty).")

  if (threshold_method == "percentile") {
    if (!is.numeric(percentile) || percentile <= 0 || percentile >= 1) {
      stop("percentile must be between 0 and 1.")
    }
    thr <- stats::quantile(mi_vals, percentile, na.rm = TRUE)
  } else {
    if (is.null(fixed_threshold) || !is.numeric(fixed_threshold) || fixed_threshold < 0) {
      stop("fixed_threshold must be a non-negative numeric value.")
    }
    thr <- fixed_threshold
  }

  adj <- matrix(0, nrow(mi_matrix), ncol(mi_matrix))
  adj[mi_matrix >= thr] <- mi_matrix[mi_matrix >= thr]
  diag(adj) <- 0
  adj <- pmax(adj, t(adj))

  rownames(adj) <- colnames(adj) <- rownames(mi_matrix)

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  ensure_weight(g)
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
#' @param percentile MI percentile threshold for graph edges
#' @return List containing graph, MI matrix, and genes
#' @export
compute_mi_network <- function(seurat_obj, cluster_id, cluster_field,
                               assay_name = "RNA",
                               counts_layer = "counts",
                               data_layer = "data",
                               top_n_genes = 5000,
                               max_cells = 1500,
                               min_expr_pct = 0.05,
                               n_cores = max(1L, parallel::detectCores() - 1L),
                               percentile = 0.95) {

  expr_subset <- get_expr_data(
    seurat_obj, cluster_id, cluster_field,
    assay_name = assay_name,
    counts_layer = counts_layer,
    data_layer = data_layer,
    top_n_genes = top_n_genes,
    max_cells = max_cells,
    min_expr_pct = min_expr_pct
  )

  mi_mat <- compute_mi_matrix_parallel(expr_subset, n_cores = n_cores)
  g <- build_mi_graph(mi_mat, threshold_method = "percentile", percentile = percentile)

  list(graph = g, mi_matrix = mi_mat, genes = rownames(expr_subset))
}

#' Percolate Graph
#'
#' Creates a consensus graph by filtering edges at a weight percentile and retaining
#' those that appear in >=freq_threshold of iterations.
#'
#' @param graph Input igraph object.
#' @param n_iterations Number of percolation iterations.
#' @param edge_fraction Fraction of edges to retain per iteration.
#' @param freq_threshold Frequency threshold for edge retention.
#' @param percentile Percentile for initial edge filtering (default: 0.95).
#' @param seed RNG seed (default: 123).
#' @return Consensus igraph object.
#' @export
percolate_graph <- function(graph,
                            n_iterations = 100,
                            edge_fraction = 0.6,
                            freq_threshold = 0.5,
                            percentile = 0.95,
                            seed = 123) {

  graph <- ensure_weight(graph)
  if (igraph::ecount(graph) == 0 || igraph::vcount(graph) == 0) return(graph)

  el_full <- igraph::as_data_frame(graph, what = "edges")
  if (!("weight" %in% colnames(el_full))) el_full$weight <- 1
  el_full <- el_full[is.finite(el_full$weight) & !is.na(el_full$weight), , drop = FALSE]
  if (nrow(el_full) == 0) return(graph)

  top_thr <- stats::quantile(el_full$weight, percentile, na.rm = TRUE)
  el_thresh <- el_full[el_full$weight >= top_thr, , drop = FALSE]

  if (nrow(el_thresh) == 0) {
    g_empty <- igraph::make_empty_graph(n = igraph::vcount(graph), directed = FALSE)
    igraph::V(g_empty)$name <- igraph::V(graph)$name
    return(g_empty)
  }

  # preserve full vertex set
  g_thresh <- igraph::graph_from_data_frame(
    el_thresh[, c("from", "to", "weight")],
    directed = FALSE,
    vertices = data.frame(name = igraph::V(graph)$name)
  )
  g_thresh <- ensure_weight(g_thresh)

  edge_list <- vector("list", n_iterations)
  for (i in seq_len(n_iterations)) {
    set.seed(seed + i)
    keep_n <- max(1L, round(igraph::ecount(g_thresh) * edge_fraction))
    keep_e <- sample(igraph::E(g_thresh), size = keep_n)
    g_i <- igraph::subgraph_from_edges(g_thresh, eids = keep_e, delete.vertices = FALSE)
    g_i <- ensure_weight(g_i)

    df <- igraph::as_data_frame(g_i, what = "edges")
    if (!("weight" %in% colnames(df))) df$weight <- 1
    df <- canonical_edge_df(df[, c("from", "to", "weight")])
    df <- df[is.finite(df$weight) & !is.na(df$weight), , drop = FALSE]
    edge_list[[i]] <- df
  }

  combined <- dplyr::bind_rows(edge_list) %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::summarise(
      freq = dplyr::n(),
      mean_weight = mean(.data$weight, na.rm = TRUE),
      .groups = "drop"
    )

  freq_thr <- ceiling(n_iterations * freq_threshold)
  cons_edges <- combined %>%
    dplyr::filter(.data$freq >= freq_thr, is.finite(.data$mean_weight))

  if (nrow(cons_edges) == 0) {
    g_empty <- igraph::make_empty_graph(n = igraph::vcount(graph), directed = FALSE)
    igraph::V(g_empty)$name <- igraph::V(graph)$name
    return(g_empty)
  }

  g_cons <- igraph::graph_from_data_frame(
    cons_edges[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = igraph::V(graph)$name)
  )

  e_df <- igraph::as_data_frame(g_cons, what = "edges")
  e_df <- canonical_edge_df(e_df[, c("from", "to")])
  e_df <- dplyr::left_join(e_df, cons_edges, by = c("from", "to"))

  w <- e_df$mean_weight
  w[!is.finite(w) | is.na(w)] <- 1
  igraph::E(g_cons)$weight <- w

  ensure_weight(g_cons)
}

#' Detect Communities
#' Detects communities using label propagation and neighborhood overlap.
#' @param graph Input igraph object.
#' @param min_size Minimum community size.
#' @param overlap_threshold Overlap threshold for adding genes (default: 0.1).
#' @param seed RNG seed.
#' @return List with groups, bridge genes, and gene2mod mapping.
#' @export
detect_communities <- function(graph, min_size = 10, overlap_threshold = 0.10, seed = 123) {

  graph <- ensure_weight(graph)

  if (igraph::ecount(graph) == 0 || igraph::vcount(graph) == 0) {
    return(list(
      groups = list(),
      bridge = character(0),
      gene2mod = tibble::tibble(Gene = character(), Module = character())
    ))
  }

  set.seed(seed)
  lp <- igraph::cluster_label_prop(graph, weights = igraph::E(graph)$weight)
  groups0 <- split(igraph::V(graph)$name, igraph::membership(lp))
  groups0 <- Filter(function(x) length(x) >= min_size, groups0)

  if (length(groups0) == 0) {
    return(list(
      groups = list(),
      bridge = character(0),
      gene2mod = tibble::tibble(Gene = character(), Module = character())
    ))
  }

  groups_exp <- groups0

  A <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
  genes <- rownames(A)
  if (is.null(genes)) genes <- igraph::V(graph)$name

  for (node in genes) {
    neighbors <- genes[A[node, ] > 0]
    if (length(neighbors) == 0) next
    for (i in seq_along(groups_exp)) {
      comm_genes <- groups_exp[[i]]
      ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
      if (ov >= overlap_threshold && !(node %in% comm_genes)) {
        groups_exp[[i]] <- c(comm_genes, node)
      }
    }
  }

  groups_exp <- lapply(groups_exp, unique)
  names(groups_exp) <- paste0("Module ", seq_along(groups_exp))

  gene2mod <- utils::stack(groups_exp)
  colnames(gene2mod) <- c("Gene", "Module")

  bridge <- dplyr::count(gene2mod, .data$Gene) %>%
    dplyr::filter(.data$n > 1) %>%
    dplyr::pull(.data$Gene)

  list(groups = groups_exp, bridge = bridge, gene2mod = gene2mod)
}

