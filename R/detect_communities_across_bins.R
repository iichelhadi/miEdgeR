#' Detect Communities Across Pseudotime Bins
#'
#' Detects communities in each pseudotime bin using percolation, Leiden clustering, and neighborhood overlap.
#'
#' @param pt_results A list of pseudotime bin results from compute_pseudotime_mi_network.
#' @param min_size Integer, minimum size of communities to retain (default: 10).
#' @param n_iterations Integer, number of percolation iterations (default: 100).
#' @param edge_fraction Numeric, fraction of edges to keep in percolation (default: 0.6).
#' @param freq_threshold Numeric, frequency threshold for consensus edges (default: 0.5).
#' @param percentile Numeric, percentile for edge thresholding (default: 0.95).
#' @param resolution_parameter Numeric, resolution parameter for Leiden clustering (default: 0.6).
#' @param overlap_threshold Numeric, threshold for neighborhood overlap expansion (default: 0.05).
#' @return A list of community detection results for each bin.
#' @export
#' @importFrom igraph cluster_leiden as_adjacency_matrix V E
detect_communities_across_bins <- function(pt_results, min_size = 10, n_iterations = 100,
                                           edge_fraction = 0.6, freq_threshold = 0.5,
                                           percentile = 0.95, resolution_parameter = 0.6,
                                           overlap_threshold = 0.05) {
  set.seed(123)  # No import needed for set.seed
  comm_results <- lapply(names(pt_results), function(bin) {
    if (!is.null(pt_results[[bin]])) {
      g_consensus <- miEdgeR::percolate_graph(pt_results[[bin]],
                                              n_iterations = n_iterations,
                                              edge_fraction = edge_fraction,
                                              freq_threshold = freq_threshold,
                                              percentile = percentile)
      lc_comm <- igraph::cluster_leiden(g_consensus, objective_function = "modularity",
                                        weights = igraph::E(g_consensus)$weight,
                                        resolution_parameter = resolution_parameter)
      lc_groups <- split(igraph::V(g_consensus)$name, igraph::membership(lc_comm))
      lc_groups <- Filter(function(x) length(x) >= min_size, lc_groups)
      adj_matrix <- igraph::as_adjacency_matrix(g_consensus, attr = "weight", sparse = FALSE)
      for (node in igraph::V(g_consensus)$name) {
        idx <- which(igraph::V(g_consensus)$name == node)
        neighbors <- igraph::V(g_consensus)$name[adj_matrix[idx, ] > 0]
        if (length(neighbors) == 0) next
        for (i in seq_along(lc_groups)) {
          comm_genes <- lc_groups[[i]]
          ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
          if (ov >= overlap_threshold && !(node %in% comm_genes)) {
            lc_groups[[i]] <- c(comm_genes, node)
          }
        }
      }
      list(bin = bin, consensus_graph = g_consensus, communities = lc_groups, sizes = sapply(lc_groups, length))
    } else {
      NULL
    }
  })
  comm_results <- Filter(Negate(is.null), comm_results)
  for (res in comm_results) {
    cat("Bin", res$bin, "- Communities (size >=", min_size, "):", length(res$communities), "\n")
    cat("Sizes:", res$sizes, "\n")
  }
  comm_results
}
