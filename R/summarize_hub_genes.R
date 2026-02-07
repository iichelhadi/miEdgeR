#' Summarize Hub Genes in Communities
#'
#' Generates a flextable summarizing the top hub genes (by degree) for each community
#' in a consensus graph.
#'
#' @param g_consensus An igraph object representing the consensus graph.
#' @param large_communities A list of communities, where each element contains the nodes of a community.
#' @param cluster_name Character, name of the cluster (default: "Cluster").
#' @param top_n Integer, number of top hub genes to include per community (default: 10).
#' @return A flextable object.
#' @export
#' @importFrom igraph induced_subgraph degree V vcount
#' @importFrom flextable flextable autofit set_header_labels bg merge_v align theme_vanilla
summarize_hub_genes <- function(g_consensus,
                               large_communities,
                               cluster_name = "Cluster",
                               top_n = 10) {

  if (!inherits(g_consensus, "igraph")) stop("g_consensus must be an igraph object.")
  if (!is.list(large_communities) || length(large_communities) == 0) {
    stop("large_communities must be a non-empty list.")
  }
  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) stop("top_n must be >= 1.")

  hub_rows <- lapply(seq_along(large_communities), function(i) {
    genes <- unique(as.character(large_communities[[i]]))
    if (length(genes) == 0) return(NULL)

    subgraph <- igraph::induced_subgraph(g_consensus, igraph::V(g_consensus)[name %in% genes])
    if (igraph::vcount(subgraph) == 0) return(NULL)

    deg <- as.numeric(unlist(igraph::degree(subgraph)))
    names(deg) <- igraph::V(subgraph)$name

    data.frame(
      Cluster   = cluster_name,
      Community = paste0("Community ", i),
      Gene      = names(deg),
      Degree    = deg,
      stringsAsFactors = FALSE
    )
  })

  hub_data <- do.call(rbind, hub_rows)

  # Handle empty result cleanly
  if (is.null(hub_data) || nrow(hub_data) == 0) {
    hub_table <- data.frame(
      Cluster = character(0),
      Community = character(0),
      Top_Hub_Genes = character(0),
      stringsAsFactors = FALSE
    )
    return(
      flextable::flextable(hub_table) |>
        flextable::autofit() |>
        flextable::set_header_labels(Cluster = "Cluster", Community = "Community", Top_Hub_Genes = "Top Hub Genes") |>
        flextable::bg(part = "header", bg = "gray") |>
        flextable::theme_vanilla()
    )
  }

  # Summarise top hubs per community (no dplyr dependency)
  hub_data <- hub_data[order(hub_data$Cluster, hub_data$Community, -hub_data$Degree), , drop = FALSE]

  split_df <- split(hub_data, paste(hub_data$Cluster, hub_data$Community, sep = "||"))
  out <- lapply(split_df, function(df) {
    df_top <- utils::head(df, top_n)
    data.frame(
      Cluster = df_top$Cluster[1],
      Community = df_top$Community[1],
      Top_Hub_Genes = paste(df_top$Gene, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  hub_table <- do.call(rbind, out)

  flextable::flextable(hub_table) |>
    flextable::autofit() |>
    flextable::set_header_labels(Cluster = "Cluster", Community = "Community", Top_Hub_Genes = "Top Hub Genes") |>
    flextable::bg(part = "header", bg = "gray") |>
    flextable::merge_v(j = "Cluster") |>
    flextable::align(j = "Cluster", align = "center", part = "body") |>
    flextable::theme_vanilla()
}

