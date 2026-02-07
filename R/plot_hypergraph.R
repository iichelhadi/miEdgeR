#' Plot Hypergraph with Overlapping Genes Highlighted
#'
#' Visualizes a hypergraph where communities are hyperedges, highlighting genes that
#' occur in >1 selected community.
#'
#' @param large_communities A list of communities, where each element contains the nodes of a community.
#' @param comm_indices Integer vector, indices of communities to plot (default: all communities).
#' @param title Character, title of the plot (default: "Hypergraph with Overlapping Communities").
#' @param key_genes Character vector, additional genes to label (default: NULL).
#' @return A ggraph plot object.
#' @export
#' @importFrom igraph graph_from_data_frame V
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 scale_size_continuous theme_void labs theme scale_fill_manual
plot_hypergraph <- function(large_communities,
                            comm_indices = NULL,
                            title = "Hypergraph with Overlapping Communities",
                            key_genes = NULL) {

  if (!is.list(large_communities) || length(large_communities) == 0) {
    stop("large_communities must be a non-empty list.")
  }

  # Default to all communities if comm_indices is NULL
  if (is.null(comm_indices)) comm_indices <- seq_along(large_communities)

  # Validate indices
  if (length(comm_indices) < 1 || any(comm_indices < 1) || any(comm_indices > length(large_communities))) {
    stop("Invalid community indices provided.")
  }

  # Subset to selected communities; coerce to unique character vectors
  selected_communities <- lapply(large_communities[comm_indices], function(x) unique(as.character(x)))

  # Overlapping genes = appear in >1 selected community (NOT intersection of all communities)
  if (length(selected_communities) < 2) {
    overlapping_genes <- character(0)
  } else {
    gene_memberships <- table(unlist(selected_communities, use.names = FALSE))
    overlapping_genes <- names(gene_memberships[gene_memberships > 1])
  }

  # Optional messaging (safe for scripts; if you don’t want console noise, remove)
  message("Overlapping genes: ", length(overlapping_genes))
  if (length(overlapping_genes) > 0) {
    message("Shared genes: ", paste(overlapping_genes, collapse = ", "))
  }

  # Build bipartite edge list: Gene <-> Community
  edges <- do.call(rbind, lapply(seq_along(selected_communities), function(i) {
    data.frame(
      Gene = selected_communities[[i]],
      Community = paste0("Comm", comm_indices[i]),
      stringsAsFactors = FALSE
    )
  }))

  g_viz <- igraph::graph_from_data_frame(edges, directed = FALSE)

  comm_nodes <- paste0("Comm", comm_indices)
  igraph::V(g_viz)$type <- igraph::V(g_viz)$name %in% comm_nodes

  # Node types
  node_type <- ifelse(
    igraph::V(g_viz)$name %in% overlapping_genes, "Overlapping",
    ifelse(igraph::V(g_viz)$type, "Community", "Gene")
  )
  igraph::V(g_viz)$node_type <- factor(node_type, levels = c("Community", "Gene", "Overlapping"))

  # Sizes
  igraph::V(g_viz)$size <- ifelse(igraph::V(g_viz)$type, 5, 2)

  # Labels: communities always; overlapping always; key_genes only if provided
  label_genes <- overlapping_genes
  if (!is.null(key_genes) && length(key_genes) > 0) {
    label_genes <- unique(c(label_genes, as.character(key_genes)))
  }
  igraph::V(g_viz)$label <- ifelse(
    igraph::V(g_viz)$type | igraph::V(g_viz)$name %in% label_genes,
    igraph::V(g_viz)$name, ""
  )

  color_values <- c("Community" = "green", "Gene" = "blue", "Overlapping" = "red")

  ggraph::ggraph(g_viz, layout = "fr") +
    ggraph::geom_edge_link(color = "grey", alpha = 0.3) +
    ggraph::geom_node_point(
      ggplot2::aes(size = .data$size, fill = .data$node_type),
      shape = 21, color = "black", show.legend = TRUE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label),
      repel = TRUE, size = 3, fontface = "bold", max.overlaps = 20
    ) +
    ggplot2::scale_size_continuous(range = c(2, 5), guide = "none") +
    ggplot2::scale_fill_manual(values = color_values, name = "Node Type", drop = TRUE) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

