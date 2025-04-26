#' Plot Enrichment Terms
#'
#' Visualizes the top enrichment terms from enrichment results in a dot plot.
#'
#' @param enrich_results A named list of enrichment results from enrich_go.
#' @param top_n_terms Integer, number of top terms to plot per community (default: 4).
#' @param title Character, title of the plot (default: "GO Terms for Communities").
#' @param ontology Character, ontology used ("BP", "CC", "MF", "ALL", "KEGG", "DO", "Reactome") (default: "BP").
#' @return A ggplot object.
#' @export
#' @importFrom dplyr bind_rows mutate
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme
plot_go_terms <- function(enrich_results, top_n_terms = 4, title = "GO Terms for Communities", ontology = "BP") {
  if (length(enrich_results) == 0) {
    warning("No enrichment results to plot.")
    return(NULL)
  }

  enrich_df <- dplyr::bind_rows(lapply(names(enrich_results), function(name) {
    head(enrich_results[[name]]@result, top_n_terms) %>% dplyr::mutate(Community = name)
  }))

  if (nrow(enrich_df) == 0) {
    warning("No enrichment terms to plot after filtering.")
    return(NULL)
  }

  ggplot(enrich_df, aes(x = Community, y = Description, size = -log10(p.adjust), color = Community)) +
    geom_point() +
    labs(x = "Community", y = paste(ontology, "Terms"), size = "-log10(p.adjust)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
