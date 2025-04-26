#' Plot Edge Gain/Loss Across Consecutive Pseudotime Bins
#'
#' Plots the count of stable, gained, and lost edges between consecutive pseudotime bins and returns edge details.
#'
#' @param pt_results A list of pseudotime bin results from compute_pseudotime_mi_network.
#' @param cluster_id Character, cluster ID to display (default: "C1").
#' @param title Character, title of the plot (default: "Edge Gain/Loss/Stable Counts Per Bin").
#' @return A list containing the ggplot object and a list of edge changes (stable, gained, lost) per bin pair.
#' @export
#' @importFrom igraph E
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal
plot_edge_gain_loss <- function(pt_results, cluster_id = "C1", title = "Edge Gain/Loss/Stable Counts Per Bin") {
  # Get bin names in order
  bin_names <- names(pt_results)
  if (length(bin_names) < 2) {
    stop("At least two bins are required to compute edge changes.")
  }

  # Compute edge changes between consecutive bins
  gl_data <- list()
  edge_details <- list()
  for (i in 1:(length(bin_names) - 1)) {
    bin1 <- bin_names[i]
    bin2 <- bin_names[i + 1]
    if (is.null(pt_results[[bin1]]) || is.null(pt_results[[bin2]])) {
      cat("Skipping", bin1, "to", bin2, ": One or both bins are NULL.\n")
      next
    }
    edge_changes <- miEdgeR::edge_gain_loss(list(pt_results[[bin1]], pt_results[[bin2]]))
    pair_name <- paste(bin1, "to", bin2)

    # Log edge counts for debugging
    cat("Edges in", bin1, ":", length(igraph::E(pt_results[[bin1]])), "\n")
    cat("Edges in", bin2, ":", length(igraph::E(pt_results[[bin2]])), "\n")
    cat("Stable edges (", pair_name, "):", length(edge_changes$stable), "\n")
    cat("Gained edges (", pair_name, "):", length(edge_changes$gained), "\n")
    cat("Lost edges (", pair_name, "):", length(edge_changes$lost), "\n")

    gl_data[[pair_name]] <- data.frame(
      Bin_Pair = pair_name,
      Cluster = cluster_id,
      Type = c("Stable", "Gained", "Lost"),
      Count = c(length(edge_changes$stable), length(edge_changes$gained), length(edge_changes$lost)),
      stringsAsFactors = FALSE
    )

    # Convert edge lists to data frames directly
    edge_details[[pair_name]] <- list(
      stable = if (length(edge_changes$stable) > 0) as.data.frame(edge_changes$stable, stringsAsFactors = FALSE) else NULL,
      gained = if (length(edge_changes$gained) > 0) as.data.frame(edge_changes$gained, stringsAsFactors = FALSE) else NULL,
      lost = if (length(edge_changes$lost) > 0) as.data.frame(edge_changes$lost, stringsAsFactors = FALSE) else NULL
    )
  }

  # Combine data for plotting
  gl_df <- do.call(rbind, gl_data)
  if (nrow(gl_df) == 0 || all(gl_df$Count == 0)) {
    warning("No edge changes detected between bins.")
    return(list(plot = NULL, edge_details = edge_details))
  }

  # Plot
  p_gl <- ggplot(gl_df, aes(x = Bin_Pair, y = Count, fill = Type)) +
    geom_col(position = "dodge") +
    labs(title = title, x = "Bin Transition", y = "Number of Edges") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(plot = p_gl, edge_details = edge_details)
}
