#' Detect Communities Across Pseudotime Bins
#'
#' Detects communities in each pseudotime bin using percolation, Leiden clustering,
#' and neighborhood overlap expansion.
#'
#' @param pt_results A named list of per-bin results from compute_pseudotime_mi_network.
#'        Each element may be an igraph object OR a list containing $graph (igraph).
#' @param min_size Integer, minimum size of communities to retain (default: 10).
#' @param n_iterations Integer, number of percolation iterations (default: 100).
#' @param edge_fraction Numeric, fraction of edges to keep in percolation (default: 0.6).
#' @param freq_threshold Numeric, frequency threshold for consensus edges (default: 0.5).
#' @param percentile Numeric, percentile for edge thresholding (default: 0.95).
#' @param resolution_parameter Numeric, resolution parameter for Leiden clustering (default: 0.6).
#' @param overlap_threshold Numeric, threshold for neighborhood overlap expansion (default: 0.05).
#' @param seed Integer RNG seed base (default: 123). Applied per bin as seed + bin_index.
#' @param verbose Logical, print progress messages (default: FALSE).
#' @return A named list with one entry per bin: consensus_graph, communities, sizes.
#' @export
#' @importFrom igraph cluster_leiden as_adjacency_matrix V E ecount vcount
detect_communities_across_bins <- function(pt_results,
                                          min_size = 10,
                                          n_iterations = 100,
                                          edge_fraction = 0.6,
                                          freq_threshold = 0.5,
                                          percentile = 0.95,
                                          resolution_parameter = 0.6,
                                          overlap_threshold = 0.05,
                                          seed = 123,
                                          verbose = FALSE) {

  stopifnot(is.list(pt_results))
  if (is.null(names(pt_results))) stop("pt_results must be a *named* list (bin names).")

  # Robust graph extractor
  .get_graph <- function(x, bin) {
    if (inherits(x, "igraph")) return(x)
    if (is.list(x) && !is.null(x$graph) && inherits(x$graph, "igraph")) return(x$graph)
    stop("pt_results[['", bin, "']] is not an igraph and does not contain $graph (igraph).")
  }
  bins <- names(pt_results)

  comm_results <- lapply(seq_along(bins), function(k) {
    bin <- bins[[k]]
    x <- pt_results[[bin]]
    if (is.null(x)) return(NULL)

    g0 <- .get_graph(x, bin)

    # IMPORTANT: inside package, call percolate_graph() directly, not miEdgeR::percolate_graph()
    set.seed(seed + k)
    g_consensus <- percolate_graph(
      g0,
      n_iterations = n_iterations,
      edge_fraction = edge_fraction,
      freq_threshold = freq_threshold,
      percentile = percentile
    )
    g_consensus <- ensure_weight(g_consensus)

    # Empty graph => no communities
    if (igraph::vcount(g_consensus) == 0 || igraph::ecount(g_consensus) == 0) {
      out <- list(
        bin = bin,
        consensus_graph = g_consensus,
        communities = list(),
        sizes = integer(0)
      )
      if (verbose) message("Bin ", bin, " - empty consensus graph (0 edges).")
      return(out)
    }

    lc_comm <- igraph::cluster_leiden(
      g_consensus,
      objective_function = "modularity",
      weights = igraph::E(g_consensus)$weight,
      resolution_parameter = resolution_parameter
    )

    lc_groups <- split(igraph::V(g_consensus)$name, igraph::membership(lc_comm))
    lc_groups <- Filter(function(v) length(v) >= min_size, lc_groups)

    if (length(lc_groups) == 0) {
      out <- list(
        bin = bin,
        consensus_graph = g_consensus,
        communities = list(),
        sizes = integer(0)
      )
      if (verbose) message("Bin ", bin, " - no communities passing min_size=", min_size, ".")
      return(out)
    }

    # Neighborhood overlap expansion
    A <- igraph::as_adjacency_matrix(g_consensus, attr = "weight", sparse = FALSE)
    genes <- rownames(A)
    if (is.null(genes)) genes <- igraph::V(g_consensus)$name

    idx_map <- seq_along(genes); names(idx_map) <- genes

    for (node in genes) {
      idx <- idx_map[[node]]
      neighbors <- genes[A[idx, ] > 0]
      if (length(neighbors) == 0) next

      for (i in seq_along(lc_groups)) {
        comm_genes <- lc_groups[[i]]
        ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
        if (ov >= overlap_threshold && !(node %in% comm_genes)) {
          lc_groups[[i]] <- c(comm_genes, node)
        }
      }
    }

    lc_groups <- lapply(lc_groups, unique)
    names(lc_groups) <- paste0("Module ", seq_along(lc_groups))
    sizes <- vapply(lc_groups, length, integer(1))

    if (verbose) {
      message("Bin ", bin, " - communities: ", length(lc_groups),
              " | sizes: ", paste(sizes, collapse = " "))
    }

    list(
      bin = bin,
      consensus_graph = g_consensus,
      communities = lc_groups,
      sizes = sizes
    )
  })

  comm_results <- Filter(Negate(is.null), comm_results)
  names(comm_results) <- vapply(comm_results, `[[`, character(1), "bin")
  comm_results
}
