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
#' @param pseudotime_field Character, metadata field containing pseudotime (default: "pseudotime").
#' @param title Character, title of the plot (default: "Module-Score Trajectories").
#' @return A ggplot object.
#' @export
#' @importFrom Seurat DefaultAssay AddModuleScore GetAssayData
#' @importFrom dplyr bind_rows
#' @importFrom ggpubr theme_pubr
#' @importFrom scales hue_pal
plot_module_trajectories <- function(comm_results,
                                     seurat_obj,
                                     bin_index,
                                     comm_indices = NULL,
                                     custom_modules = NULL,
                                     cluster_field = "RNA_snn_res.0.1",
                                     cluster_id = "1",
                                     assay = "RNA",
                                     pseudotime_field = "pseudotime",
                                     title = "Module-Score Trajectories") {

  if (!is.list(comm_results) || length(comm_results) == 0) stop("comm_results must be a non-empty list.")
  if (!bin_index %in% seq_along(comm_results)) stop("bin_index must be a valid index within comm_results.")

  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.")
  if (!cluster_field %in% colnames(seurat_obj@meta.data)) stop("cluster_field not found in seurat_obj@meta.data.")
  if (!pseudotime_field %in% colnames(seurat_obj@meta.data)) stop("pseudotime_field not found in seurat_obj@meta.data.")

  # Extract communities for the specified bin
  if (is.null(comm_results[[bin_index]]$communities) || !is.list(comm_results[[bin_index]]$communities)) {
    stop("comm_results[[bin_index]]$communities must be a list.")
  }
  large_communities <- lapply(comm_results[[bin_index]]$communities, function(x) unique(as.character(unlist(x))))

  # Build module list
  if (!is.null(custom_modules)) {
    if (!is.list(custom_modules) || is.null(names(custom_modules))) stop("custom_modules must be a named list.")
    mods <- lapply(custom_modules, function(x) unique(as.character(x)))
  } else {
    if (is.null(comm_indices)) comm_indices <- seq_along(large_communities)
    if (length(comm_indices) < 1 || any(comm_indices < 1) || any(comm_indices > length(large_communities))) {
      stop("Invalid community indices provided.")
    }
    mods <- setNames(
      lapply(comm_indices, function(i) large_communities[[i]]),
      paste0(as.character(cluster_id), "_M", comm_indices)
    )
  }

  # Prepare object + filter cells once
  so <- seurat_obj
  Seurat::DefaultAssay(so) <- assay

  cluster_vec <- as.character(so@meta.data[[cluster_field]])
  pt_vec <- so@meta.data[[pseudotime_field]]
  keep_cells <- colnames(so)[cluster_vec == as.character(cluster_id) & !is.na(pt_vec)]

  if (length(keep_cells) == 0) stop("No cells remain after filtering by cluster_field/cluster_id and non-NA pseudotime.")

  # Ensure features exist in assay
  feat_universe <- rownames(Seurat::GetAssayData(so, assay = assay, slot = "data"))

  # Add module scores (skip empty modules)
  for (modName in names(mods)) {
    feats <- intersect(mods[[modName]], feat_universe)
    if (length(feats) == 0) {
      warning("Skipping module '", modName, "': 0 genes found in assay.")
      next
    }
    so <- Seurat::AddModuleScore(so, features = list(feats), name = modName)
  }

  valid_mods <- names(mods)[paste0(names(mods), "1") %in% colnames(so@meta.data)]
  if (length(valid_mods) == 0) stop("No module score columns were created (all modules empty or missing).")

  # Trajectory dataframe
  traj_df <- dplyr::bind_rows(lapply(valid_mods, function(modName) {
    score_col <- paste0(modName, "1")
    data.frame(
      Cell       = keep_cells,
      Cluster    = as.character(cluster_id),
      Module     = modName,
      Pseudotime = so@meta.data[keep_cells, pseudotime_field],
      Score      = so@meta.data[keep_cells, score_col],
      stringsAsFactors = FALSE
    )
  }))

  # Labels + colors
  module_labels <- stats::setNames(valid_mods, valid_mods)
  module_color <- scales::hue_pal()(length(valid_mods))
  names(module_color) <- valid_mods

  ggplot2::ggplot(traj_df, ggplot2::aes(x = .data$Pseudotime, y = .data$Score, color = .data$Module, fill = .data$Module)) +
    ggplot2::geom_smooth(se = TRUE, span = 0.5) +
    ggplot2::facet_wrap(~ Cluster, scales = "free_x") +
    ggplot2::labs(title = title, x = "Pseudotime", y = "Module Score") +
    ggpubr::theme_pubr(legend = "right") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 18)
    ) +
    ggplot2::scale_color_manual(values = module_color, labels = module_labels) +
    ggplot2::scale_fill_manual(values = module_color, labels = module_labels)
}

