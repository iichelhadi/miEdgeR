#' Plot MI Network with Community Coloring
#'
#' Visualizes a subset of an MI network with community coloring, labeling top nodes from each community.
#'
#' @param g_consensus An igraph object representing the consensus graph.
#' @param large_communities A list of communities, where each element contains the nodes of a community.
#' @param sample_size Integer, number of nodes to sample from the graph (default: 100).
#' @param top_n_per_comm Integer, number of top nodes per community to label (default: 5).
#' @param title Character, title of the plot (default: "Network Communities").
#' @return A ggraph plot object.
#' @export
#' @importFrom igraph induced_subgraph degree V
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_alpha
#' @importFrom ggplot2 scale_size_continuous theme_void labs theme scale_fill_manual
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
plot_mi_network <- function(g_consensus,
                            large_communities,
                            sample_size = 100,
                            top_n_per_comm = 5,
                            title = "Network Communities") {

  if (!inherits(g_consensus, "igraph")) stop("g_consensus must be an igraph object.")
  if (!is.list(large_communities) || length(large_communities) == 0) stop("large_communities must be a non-empty list.")

  # Ensure vertex names exist
  if (is.null(igraph::V(g_consensus)$name)) {
    igraph::V(g_consensus)$name <- as.character(seq_len(igraph::vcount(g_consensus)))
  }

  # Ensure edge weights exist
  if (!("weight" %in% igraph::edge_attr_names(g_consensus))) {
    igraph::E(g_consensus)$weight <- 1
  } else {
    w <- igraph::E(g_consensus)$weight
    bad <- !is.finite(w) | is.na(w)
    if (any(bad)) igraph::E(g_consensus)$weight[bad] <- 1
  }

  set.seed(123)
  n <- igraph::vcount(g_consensus)
  sample_size <- max(2, min(as.integer(sample_size), n))  # need >=2 nodes for a sensible plot

  sampled_nodes <- sample(igraph::V(g_consensus)$name, sample_size)
  subgraph <- igraph::induced_subgraph(g_consensus, sampled_nodes)

  # Node size = degree in subgraph
  igraph::V(subgraph)$size <- igraph::degree(subgraph)

  # Assign community membership (by gene name)
  comm_membership <- rep(NA_integer_, length(igraph::V(g_consensus)$name))
  names(comm_membership) <- igraph::V(g_consensus)$name

  for (i in seq_along(large_communities)) {
    comm_membership[as.character(large_communities[[i]])] <- i
  }

  igraph::V(subgraph)$community <- comm_membership[igraph::V(subgraph)$name]

  # Select top nodes per community by degree
  deg_filtered <- igraph::degree(subgraph)
  comm_ids <- sort(unique(stats::na.omit(igraph::V(subgraph)$community)))

  top_nodes <- unlist(lapply(comm_ids, function(comm) {
    comm_nodes <- igraph::V(subgraph)$name[igraph::V(subgraph)$community == comm]
    if (length(comm_nodes) == 0) return(character(0))
    deg_subset <- deg_filtered[comm_nodes]
    names(sort(deg_subset, decreasing = TRUE))[seq_len(min(top_n_per_comm, length(deg_subset)))]
  }), use.names = FALSE)

  igraph::V(subgraph)$label <- ifelse(igraph::V(subgraph)$name %in% top_nodes, igraph::V(subgraph)$name, NA)

  # Color palette for communities (only those that exist in large_communities list)
  n_communities <- length(large_communities)
  if (n_communities <= 2) {
    community_colors <- c("#f56642", "#00B6EB")[seq_len(n_communities)]
  } else {
    base_cols <- RColorBrewer::brewer.pal(min(n_communities, 8), "Set2")
    community_colors <- if (n_communities > 8) grDevices::colorRampPalette(base_cols)(n_communities) else base_cols
  }
  names(community_colors) <- as.character(seq_len(n_communities))
  community_labels <- paste("Community", seq_len(n_communities))

  ggraph::ggraph(subgraph, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_alpha = .data$weight, edge_width = .data$weight),
      color = "grey60", show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = .data$size, fill = factor(.data$community)),
      shape = 21, color = "black", show.legend = TRUE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label),
      fontface = "bold", repel = TRUE, size = 6, color = "black"
    ) +
    ggplot2::scale_size_continuous(range = c(3, 15), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.2, 1)) +
    ggplot2::scale_fill_manual(values = community_colors, labels = community_labels, name = "Community") +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

