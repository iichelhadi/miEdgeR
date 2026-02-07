#' Summarize Hub Genes Across Pseudotime Bins
#'
#' Generates a flextable summarizing the top hub genes (by degree) for each community
#' across pseudotime bins.
#'
#' @param comm_results A list of community detection results from detect_communities_across_bins.
#' @param pt_results A list of pseudotime bin graphs from compute_pseudotime_mi_network.
#' @param top_n Integer, number of top hub genes to include per community (default: 10).
#' @return A flextable object.
#' @export
#' @importFrom igraph induced_subgraph degree V vcount
#' @importFrom flextable flextable autofit set_header_labels bg merge_v align theme_vanilla
summarize_hub_genes_across_bins <- function(comm_results, pt_results, top_n = 10) {

  if (!is.list(comm_results) || length(comm_results) == 0) {
    stop("comm_results must be a non-empty list.")
  }
  if (!is.list(pt_results) || length(pt_results) == 0) {
    stop("pt_results must be a non-empty list.")
  }
  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
    stop("top_n must be a single integer >= 1.")
  }

  hub_tables <- lapply(comm_results, function(res) {
    bin <- res$bin
    g_bin <- pt_results[[bin]]

    if (is.null(g_bin) || !inherits(g_bin, "igraph")) return(NULL)
    if (is.null(res$communities) || length(res$communities) == 0) return(NULL)

    rows <- lapply(seq_along(res$communities), function(i) {
      genes <- unique(as.character(res$communities[[i]]))
      if (length(genes) == 0) return(NULL)

      subgraph <- igraph::induced_subgraph(g_bin, igraph::V(g_bin)[name %in% genes])
      if (igraph::vcount(subgraph) == 0) return(NULL)

      deg <- as.numeric(unlist(igraph::degree(subgraph)))
      names(deg) <- igraph::V(subgraph)$name

      data.frame(
        Bin = bin,
        Community = paste0("Community ", i),
        Gene = names(deg),
        Degree = deg,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  })

  hub_tables <- do.call(rbind, hub_tables)

  # Empty-safe return
  if (is.null(hub_tables) || nrow(hub_tables) == 0) {
    out <- data.frame(
      Bin = character(0),
      Community = character(0),
      Top_Hub_Genes = character(0),
      stringsAsFactors = FALSE
    )
    return(
      flextable::flextable(out) |>
        flextable::autofit() |>
        flextable::set_header_labels(Bin = "Bin", Community = "Community", Top_Hub_Genes = "Top Hub Genes") |>
        flextable::bg(part = "header", bg = "gray") |>
        flextable::theme_vanilla()
    )
  }

  # Base-R summarize top hubs per (Bin, Community)
  hub_tables <- hub_tables[order(hub_tables$Bin, hub_tables$Community, -hub_tables$Degree), , drop = FALSE]
  key <- paste(hub_tables$Bin, hub_tables$Community, sep = "||")
  split_df <- split(hub_tables, key)

  out <- lapply(split_df, function(df) {
    df_top <- utils::head(df, top_n)
    data.frame(
      Bin = df_top$Bin[1],
      Community = df_top$Community[1],
      Top_Hub_Genes = paste(df_top$Gene, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  hub_summary <- do.call(rbind, out)

  flextable::flextable(hub_summary) |>
    flextable::autofit() |>
    flextable::set_header_labels(Bin = "Bin", Community = "Community", Top_Hub_Genes = "Top Hub Genes") |>
    flextable::bg(part = "header", bg = "gray") |>
    flextable::merge_v(j = "Bin") |>
    flextable::align(j = "Bin", align = "center", part = "body") |>
    flextable::theme_vanilla()
}

