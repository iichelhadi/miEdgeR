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
plot_mi_network <- function(g_consensus, large_communities, sample_size = 100, top_n_per_comm = 5, title = "Network Communities") {
  set.seed(123)
  subgraph <- igraph::induced_subgraph(g_consensus, sample(igraph::V(g_consensus)$name, sample_size))
  igraph::V(subgraph)$size <- igraph::degree(subgraph)

  # Assign community membership
  comm_membership <- rep(NA, length(igraph::V(g_consensus)$name))
  names(comm_membership) <- igraph::V(g_consensus)$name
  for (i in seq_along(large_communities)) {
    comm_membership[large_communities[[i]]] <- i
  }
  igraph::V(subgraph)$community <- comm_membership[igraph::V(subgraph)$name]

  # Select top nodes per community by degree
  deg_filtered <- igraph::degree(subgraph)
  top_nodes <- unlist(lapply(unique(na.omit(igraph::V(subgraph)$community)), function(comm) {
    comm_nodes <- igraph::V(subgraph)$name[igraph::V(subgraph)$community == comm]
    if (length(comm_nodes) == 0) return(NULL)
    deg_subset <- deg_filtered[comm_nodes]
    names(sort(deg_subset, decreasing = TRUE))[1:min(top_n_per_comm, length(deg_subset))]
  }))
  igraph::V(subgraph)$label <- ifelse(igraph::V(subgraph)$name %in% top_nodes, igraph::V(subgraph)$name, NA)

  # Dynamic color palette for communities
  n_communities <- length(large_communities)
  if (n_communities <= 2) {
    community_colors <- c("#f56642", "#00B6EB")[1:n_communities]
  } else {
    community_colors <- RColorBrewer::brewer.pal(min(n_communities, 8), "Set2")
    if (n_communities > 8) {
      community_colors <- grDevices::colorRampPalette(community_colors)(n_communities)
    }
  }
  names(community_colors) <- as.character(1:n_communities)
  community_labels <- paste("Community", 1:n_communities)

  ggraph::ggraph(subgraph, layout = "fr") +
    ggraph::geom_edge_link(aes(edge_alpha = weight, edge_width = weight), color = "grey60", show.legend = FALSE) +
    ggraph::geom_node_point(aes(size = size, fill = factor(community)), shape = 21, color = "black", show.legend = TRUE) +
    ggraph::geom_node_text(aes(label = label), fontface = "bold", repel = TRUE, size = 6, color = "black") +
    ggplot2::scale_size_continuous(range = c(3, 15), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.2, 1)) +
    ggplot2::scale_fill_manual(values = community_colors, labels = community_labels, name = "Community") +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
}
