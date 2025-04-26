#' Plot Hypergraph with Overlapping Genes Highlighted
#'
#' Visualizes a hypergraph where communities are hyperedges, highlighting overlapping genes between specified communities.
#'
#' @param large_communities A list of communities, where each element contains the nodes of a community.
#' @param comm_indices Integer vector, indices of communities to plot (default: all communities).
#' @param title Character, title of the plot (default: "Hypergraph with Overlapping Communities").
#' @param key_genes Character vector, additional genes to label (default: c("PRSS1", "PRSS2")).
#' @return A ggraph plot object.
#' @export
#' @importFrom igraph graph_from_data_frame V
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 scale_size_continuous theme_void labs theme scale_fill_manual
plot_hypergraph <- function(large_communities, comm_indices = NULL, title = "Hypergraph with Overlapping Communities", key_genes = c("PRSS1", "PRSS2")) {
  # Default to all communities if comm_indices is NULL
  if (is.null(comm_indices)) {
    comm_indices <- seq_along(large_communities)
  }

  # Validate indices
  if (length(comm_indices) < 1 || any(comm_indices > length(large_communities))) {
    stop("Invalid community indices provided.")
  }

  # Subset to selected communities
  selected_communities <- large_communities[comm_indices]

  # Identify overlapping genes across selected communities
  if (length(selected_communities) < 2) {
    overlapping_genes <- character(0)
  } else {
    overlapping_genes <- Reduce(intersect, selected_communities)
  }
  cat("Overlapping genes:", length(overlapping_genes), "\n")
  if (length(overlapping_genes) > 0) {
    cat("Shared genes:", paste(overlapping_genes, collapse = ", "), "\n")
  }

  # Create hypergraph with selected communities
  groups <- selected_communities
  edges <- data.frame()
  for (i in seq_along(groups)) {
    edges <- rbind(edges, data.frame(Gene = groups[[i]], Community = paste("Comm", comm_indices[i], sep = "")))
  }
  g_viz <- igraph::graph_from_data_frame(edges, directed = FALSE)
  V(g_viz)$type <- V(g_viz)$name %in% paste0("Comm", comm_indices)
  V(g_viz)$color <- ifelse(V(g_viz)$name %in% overlapping_genes, "Overlapping", ifelse(V(g_viz)$type, "Community", "Gene"))
  V(g_viz)$size <- ifelse(V(g_viz)$type, 5, 2)
  V(g_viz)$label <- ifelse(V(g_viz)$type | V(g_viz)$name %in% c(overlapping_genes, key_genes), V(g_viz)$name, "")

  # Dynamic color palette for nodes
  unique_colors <- unique(V(g_viz)$color)
  color_values <- c("Community" = "green", "Gene" = "blue", "Overlapping" = "red")[unique_colors]

  ggraph(g_viz, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.3) +
    geom_node_point(aes(size = size, fill = color), shape = 21, color = "black", show.legend = TRUE) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold", max.overlaps = 20) +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    scale_fill_manual(values = color_values, name = "Node Type") +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
}
