#' Summarize Hub Genes in Communities
#'
#' Generates a flextable summarizing the top hub genes (by degree) for each community in a consensus graph.
#'
#' @param g_consensus An igraph object representing the consensus graph.
#' @param large_communities A list of communities, where each element contains the nodes of a community.
#' @param cluster_name Character, name of the cluster (default: "Cluster").
#' @param top_n Integer, number of top hub genes to include per community (default: 10).
#' @return A flextable object.
#' @export
#' @importFrom igraph induced_subgraph degree V
#' @importFrom dplyr group_by arrange desc slice_head summarise
#' @importFrom flextable flextable autofit set_header_labels bg merge_v align theme_vanilla
summarize_hub_genes <- function(g_consensus, large_communities, cluster_name = "Cluster", top_n = 10) {
  hub_data <- lapply(seq_along(large_communities), function(i) {
    genes <- large_communities[[i]]
    subgraph <- igraph::induced_subgraph(g_consensus, igraph::V(g_consensus)[name %in% genes])
    if (igraph::vcount(subgraph) == 0) return(NULL)
    deg <- igraph::degree(subgraph)
    data.frame(
      Cluster = cluster_name,
      Community = paste("Community", i),
      Gene = names(deg),
      Degree = deg,
      stringsAsFactors = FALSE
    )
  })

  hub_data <- do.call(rbind, hub_data)
  hub_table <- hub_data %>%
    dplyr::group_by(Cluster, Community) %>%
    dplyr::arrange(dplyr::desc(Degree)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::summarise(Top_Hub_Genes = paste(Gene, collapse = ", "), .groups = "drop")

  ft <- flextable::flextable(hub_table) %>%
    flextable::autofit() %>%
    flextable::set_header_labels(Cluster = "Cluster", Community = "Community", Top_Hub_Genes = "Top Hub Genes") %>%
    flextable::bg(part = "header", bg = "gray") %>%
    flextable::merge_v(j = "Cluster") %>%
    flextable::align(j = "Cluster", align = "center", part = "body") %>%
    flextable::theme_vanilla()
  ft
}
