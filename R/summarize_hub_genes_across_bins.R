#' Summarize Hub Genes Across Pseudotime Bins
#'
#' Generates a flextable summarizing the top hub genes (by degree) for each community across pseudotime bins.
#'
#' @param comm_results A list of community detection results from detect_communities_across_bins.
#' @param pt_results A list of pseudotime bin results from compute_pseudotime_mi_network.
#' @param top_n Integer, number of top hub genes to include per community (default: 10).
#' @return A flextable object.
#' @export
#' @importFrom igraph induced_subgraph degree V
#' @importFrom dplyr group_by arrange desc slice_head summarise
#' @importFrom flextable flextable autofit set_header_labels bg merge_v align theme_vanilla
summarize_hub_genes_across_bins <- function(comm_results, pt_results, top_n = 10) {
  hub_tables <- lapply(comm_results, function(res) {
    bin <- res$bin
    hub_data <- lapply(seq_along(res$communities), function(i) {
      genes <- res$communities[[i]]
      subgraph <- igraph::induced_subgraph(pt_results[[bin]], igraph::V(pt_results[[bin]])[name %in% genes])
      if (igraph::vcount(subgraph) == 0) return(NULL)
      deg <- igraph::degree(subgraph)
      data.frame(
        Bin = bin,
        Group = paste("Group", i),
        Gene = names(deg),
        Degree = deg,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, hub_data)
  })
  hub_tables <- do.call(rbind, hub_tables)
  hub_summary <- hub_tables %>%
    dplyr::group_by(Bin, Group) %>%
    dplyr::arrange(dplyr::desc(Degree)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::summarise(Top_Hub_Genes = paste(Gene, collapse = ", "), .groups = "drop")

  ft <- flextable::flextable(hub_summary) %>%
    flextable::autofit() %>%
    flextable::set_header_labels(Bin = "Bin", Group = "Group", Top_Hub_Genes = "Top Hub Genes") %>%
    flextable::bg(part = "header", bg = "gray") %>%
    flextable::merge_v(j = "Bin") %>%
    flextable::align(j = "Bin", align = "center", part = "body") %>%
    flextable::theme_vanilla()
  ft
}
