#' Plot Module-Score Trajectories Across Pseudotime
#'
#' Computes and plots module scores for specified communities and genes along pseudotime.
#'
#' @param comm_results A list of community detection results from detect_communities_across_bins.
#' @param seurat_obj A Seurat object containing pseudotime and cluster metadata.
#' @param bin_index Integer, index of the pseudotime bin to use (e.g., 2 for PT2).
#' @param comm_indices Integer vector, indices of communities to plot (default: NULL, uses all communities).
#' @param custom_modules A named list of custom gene sets to plot (default: NULL; if provided, overrides comm_indices).
#' @param cluster_field Character, metadata field for cluster IDs (default: "RNA_snn_res.0.1").
#' @param cluster_id Character, cluster ID to filter cells (default: "1").
#' @param assay Character, assay to use for module scoring (default: "RNA").
#' @param title Character, title of the plot (default: "Module-Score Trajectories").
#' @return A ggplot object.
#' @export
#' @importFrom Seurat DefaultAssay AddModuleScore GetAssayData
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_smooth facet_wrap labs geom_hline theme scale_color_manual scale_fill_manual
#' @importFrom ggpubr theme_pubr
#' @importFrom scales hue_pal
plot_module_trajectories <- function(comm_results, seurat_obj, bin_index, comm_indices = NULL, custom_modules = NULL,
                                     cluster_field = "RNA_snn_res.0.1", cluster_id = "1", assay = "RNA",
                                     title = "Module-Score Trajectories") {
  # Validate inputs
  if (!bin_index %in% seq_along(comm_results)) {
    stop("bin_index must be a valid index within comm_results.")
  }

  # Extract communities for the specified bin
  large_communities <- lapply(comm_results[[bin_index]]$communities, unlist)

  # If custom_modules is provided, use it exclusively; otherwise, use comm_indices
  if (!is.null(custom_modules)) {
    if (!is.list(custom_modules) || is.null(names(custom_modules))) {
      stop("custom_modules must be a named list.")
    }
    mods <- custom_modules
  } else {
    # Default to all communities if comm_indices is NULL
    if (is.null(comm_indices)) {
      comm_indices <- seq_along(large_communities)
    }

    # Validate comm_indices
    if (length(comm_indices) < 1 || any(comm_indices > length(large_communities))) {
      stop("Invalid community indices provided.")
    }

    # Create modules for selected communities with dynamic naming
    mods <- lapply(comm_indices, function(i) {
      setNames(list(large_communities[[i]]), paste0(as.character(cluster_id), "_M", i))
    })
    mods <- do.call(c, mods)
  }

  # Compute module scores
  so <- seurat_obj
  DefaultAssay(so) <- assay
  for (modName in names(mods)) {
    so <- AddModuleScore(so, features = list(mods[[modName]]), name = modName)
  }

  # Create trajectory data frame
  traj_df <- dplyr::bind_rows(lapply(names(mods), function(modName) {
    keep_cells <- colnames(so)[so@meta.data[[cluster_field]] == cluster_id & !is.na(so@meta.data$pseudotime)]
    score_col <- paste0(modName, "1")
    data.frame(
      Cell = keep_cells,
      Cluster = cluster_id,
      Module = modName,
      Pseudotime = so@meta.data[keep_cells, "pseudotime"],
      Score = so@meta.data[keep_cells, score_col],
      stringsAsFactors = FALSE
    )
  }))

  # Define module labels and colors
  n_modules <- length(mods)
  module_labels <- setNames(names(mods), names(mods))
  module_color <- scales::hue_pal()(n_modules)
  names(module_color) <- names(mods)

  # Plot trajectories
  p_modtraj <- ggplot(traj_df, aes(x = Pseudotime, y = Score, color = Module, fill = Module)) +
    geom_smooth(se = TRUE, span = 0.5) +
    facet_wrap(~ Cluster, scales = "free_x") +
    labs(title = title, x = "Pseudotime", y = "Module Score") +
    ggpubr::theme_pubr(legend = "right") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 18)) +
    scale_color_manual(values = module_color, labels = module_labels) +
    scale_fill_manual(values = module_color, labels = module_labels)

  p_modtraj
}
