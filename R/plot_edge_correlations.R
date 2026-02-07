#' Plot Edge-Pseudotime Correlation Density
#'
#' Computes and plots the density of edge–pseudotime correlations for a specified pseudotime bin.
#'
#' @param comm_results Output from detect_communities_across_bins()
#' @param pt_results Output from compute_pseudotime_mi_network()
#' @param seurat_obj Seurat object with pseudotime metadata
#' @param bin_name Pseudotime bin name (e.g. "PT2")
#' @param cluster_field Metadata column for cluster labels
#' @param cluster_id Cluster ID to subset
#' @param assay Assay name
#' @param layer Data layer
#' @param pseudotime_col Pseudotime column name
#' @param method Correlation method
#' @param edge_threshold MI edge threshold
#' @param n_cores Number of cores
#' @param include_cross Logical; include cross-community edges (default FALSE)
#' @param title Optional plot title
#'
#' @return ggplot object or NULL
#' @export
plot_edge_correlations <- function(
    comm_results, pt_results, seurat_obj, bin_name,
    cluster_field = "RNA_snn_res.0.1", cluster_id = "1",
    assay = "RNA", layer = "data", pseudotime_col = "pseudotime",
    method = "spearman", edge_threshold = 0.95, n_cores = 5,
    include_cross = FALSE,
    title = NULL
) {

  stopifnot(bin_name %in% names(pt_results))

  ## ---- get graph ----
  pt_obj <- pt_results[[bin_name]]
  g_bin <- if (inherits(pt_obj, "igraph")) pt_obj else pt_obj$graph
  stopifnot(inherits(g_bin, "igraph"))

  ## ---- subset cells ----
  keep_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_field]] == cluster_id]
  stopifnot(length(keep_cells) > 0)
  so_sub <- seurat_obj[, keep_cells]

  ## ---- edge–pseudotime correlation ----
  ec_df <- edge_pt_correlation(
    pt_obj,
    so_sub,
    pseudotime_col = pseudotime_col,
    assay = assay,
    layer = layer,
    method = method,
    edge_threshold = edge_threshold,
    n_cores = n_cores
  )

  if (is.null(ec_df) || nrow(ec_df) == 0) return(NULL)

  ec_df <- ec_df[is.finite(ec_df$cor_pt), , drop = FALSE]
  if (nrow(ec_df) == 0) return(NULL)

  ## ---- communities for this bin ----
  comm_map <- setNames(comm_results, vapply(comm_results, `[[`, "", "bin"))
  comms <- comm_map[[bin_name]]$communities
  stopifnot(is.list(comms), length(comms) > 0)

  vertex_names <- igraph::V(g_bin)$name

  gene2comm <- do.call(rbind, lapply(seq_along(comms), function(i) {
    data.frame(
      gene = intersect(comms[[i]], vertex_names),
      community = paste0("Community ", i),
      stringsAsFactors = FALSE
    )
  }))

  ## ---- annotate edges ----
  ec_df <- merge(ec_df, gene2comm, by.x = "from", by.y = "gene")
  ec_df <- merge(ec_df, gene2comm, by.x = "to", by.y = "gene",
                 suffixes = c("_from", "_to"))

  ec_df$Community <- ifelse(
    ec_df$community_from == ec_df$community_to,
    ec_df$community_from,
    "Cross-Community"
  )

  if (!include_cross) {
    ec_df <- ec_df[ec_df$Community != "Cross-Community", , drop = FALSE]
  }

  if (nrow(ec_df) == 0) return(NULL)

  ec_df$Community <- factor(ec_df$Community)

  ## ---- colors ----
  comm_levels <- levels(ec_df$Community)
  base_cols <- c("Community 1" = "#f56642", "Community 2" = "#00B6EB")
  color_values <- base_cols[comm_levels]
  color_values[is.na(color_values)] <- "grey"

  if (is.null(title)) {
    title <- paste0(
      "Edge–Pseudotime Correlation Density (",
      bin_name, ", edge_threshold=", edge_threshold, ")"
    )
  }

  ## ---- plot ----
  ggplot2::ggplot(ec_df, ggplot2::aes(x = cor_pt, color = Community, fill = Community)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
    ggplot2::scale_fill_manual(values = color_values, drop = TRUE) +
    ggplot2::labs(title = title, x = "Spearman ρ", y = "Density") +
    ggplot2::theme_minimal()
}
