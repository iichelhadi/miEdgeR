#' Plot Edge-Pseudotime Correlation Density
#'
#' Computes and plots the density of edge-pseudotime correlations for a specified pseudotime bin and cluster.
#'
#' @param comm_results A list of community detection results from detect_communities_across_bins.
#' @param pt_results A list of pseudotime bin results from compute_pseudotime_mi_network.
#' @param seurat_obj A Seurat object containing pseudotime and cluster metadata.
#' @param bin_name Character, name of the pseudotime bin to use (e.g., "PT2").
#' @param cluster_field Character, metadata field for cluster IDs (default: "RNA_snn_res.0.1").
#' @param cluster_id Character, cluster ID to filter cells (default: "1").
#' @param assay Character, assay to use (default: "RNA").
#' @param layer Character, data layer to use (default: "data").
#' @param pseudotime_col Character, metadata column for pseudotime (default: "pseudotime").
#' @param method Character, correlation method (default: "spearman").
#' @param edge_threshold Numeric, edge weight threshold (default: 0.95).
#' @param n_cores Integer, number of cores for parallel computation (default: 5).
#' @param title Character, title of the plot (default: dynamically generated).
#' @return A ggplot object.
#' @export
#' @importFrom igraph V
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_density scale_color_manual scale_fill_manual labs theme_minimal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
plot_edge_correlations <- function(comm_results, pt_results, seurat_obj, bin_name,
                                   cluster_field = "RNA_snn_res.0.1", cluster_id = "1",
                                   assay = "RNA", layer = "data", pseudotime_col = "pseudotime",
                                   method = "spearman", edge_threshold = 0.95, n_cores = 5,
                                   title = NULL) {
  # Assign names to comm_results entries
  names(comm_results) <- sapply(comm_results, function(res) res$bin)

  # Check if comm_results contains the selected bin
  if (length(comm_results) == 0 || !bin_name %in% names(comm_results)) {
    stop("No valid communities found for bin ", bin_name, ". Check pseudotime binning or community detection.")
  }

  # Compute edge-pseudotime correlations
  ec_df <- edge_pt_correlation(pt_results[[bin_name]],
                               seurat_obj[, seurat_obj@meta.data[[cluster_field]] == cluster_id],
                               pseudotime_col = pseudotime_col, assay = assay, layer = layer,
                               method = method, edge_threshold = edge_threshold, n_cores = n_cores)

  # Check for non-finite values in cor_pt and filter them out
  if (any(!is.finite(ec_df$cor_pt))) {
    warning("Non-finite values found in cor_pt. Filtering them out.")
    ec_df <- ec_df[is.finite(ec_df$cor_pt), ]
  }

  # Assign communities to edges dynamically
  large_communities <- lapply(comm_results[[which(names(comm_results) == bin_name)]]$communities, unlist)
  comm_membership <- rep(NA, length(igraph::V(pt_results[[bin_name]])$name))
  names(comm_membership) <- igraph::V(pt_results[[bin_name]])$name
  for (i in seq_along(large_communities)) {
    comm_membership[large_communities[[i]]] <- paste("Community", i)
  }

  ec_df$Community <- sapply(1:nrow(ec_df), function(i) {
    from_comm <- comm_membership[ec_df$from[i]]
    to_comm <- comm_membership[ec_df$to[i]]
    if (is.na(from_comm) || is.na(to_comm)) return(NA)
    if (from_comm == to_comm) return(from_comm)
    return("Cross-Community")
  })

  # Filter out NA values in Community
  ec_df <- ec_df[!is.na(ec_df$Community), ]
  if (nrow(ec_df) == 0) {
    warning("No valid data to plot after filtering.")
    return(NULL)
  }
  ec_df$Community <- factor(ec_df$Community, levels = unique(ec_df$Community))

  # Dynamic color palette for communities
  n_communities <- sum(grepl("^Community", levels(ec_df$Community)))
  if (n_communities <= 2) {
    community_colors <- c("Community 1" = "#f56642", "Community 2" = "#00B6EB")
  } else {
    community_colors <- RColorBrewer::brewer.pal(n_communities, "Set2")
    if (n_communities > 8) {
      community_colors <- grDevices::colorRampPalette(community_colors)(n_communities)
    }
    names(community_colors) <- paste("Community", 1:n_communities)
  }
  color_values <- c(community_colors, "Cross-Community" = "grey")[levels(ec_df$Community)]

  # Set default title if not provided
  if (is.null(title)) {
    title <- paste("Edge-Pseudotime Correlation Density (", bin_name, ", Top 5% Edges)")
  }

  # Plot density per community
  p_ec <- ggplot2::ggplot(ec_df, aes(x = cor_pt, color = Community, fill = Community)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
    ggplot2::scale_fill_manual(values = color_values, drop = TRUE) +
    ggplot2::labs(title = title, x = "Spearman ρ", y = "Density") +
    ggplot2::theme_minimal()

  p_ec
}
