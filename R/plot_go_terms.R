#' Plot Enrichment Terms
#'
#' Visualizes the top enrichment terms from enrichment results in a dot plot.
#'
#' @param enrich_results A named list of enrichment results from enrich_go (clusterProfiler enrichResult
#'   objects OR data.frames with columns Description and p.adjust).
#' @param top_n_terms Integer, number of top terms to plot per community (default: 4).
#' @param title Character, title of the plot (default: "GO Terms for Communities").
#' @param ontology Character, ontology used ("BP", "CC", "MF", "ALL", "KEGG", "DO", "Reactome") (default: "BP").
#' @return A ggplot object or NULL.
#' @export
plot_go_terms <- function(enrich_results,
                          top_n_terms = 4,
                          title = "GO Terms for Communities",
                          ontology = "BP") {

  if (is.null(enrich_results) || length(enrich_results) == 0) {
    warning("No enrichment results to plot.")
    return(NULL)
  }

  enrich_df <- dplyr::bind_rows(lapply(names(enrich_results), function(nm) {
    er <- enrich_results[[nm]]
    df <- NULL

    # clusterProfiler enrichResult / compareClusterResult typically has @result
    if (!is.null(er@result)) {
      df <- er@result
    } else if (is.data.frame(er)) {
      df <- er
    } else {
      df <- tryCatch(as.data.frame(er), error = function(e) NULL)
    }

    if (is.null(df) || nrow(df) == 0) return(NULL)

    # minimal required columns
    if (!all(c("Description", "p.adjust") %in% colnames(df))) return(NULL)

    df <- df[seq_len(min(top_n_terms, nrow(df))), , drop = FALSE]
    df$Community <- nm
    df
  }))

  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    warning("No enrichment terms to plot after filtering.")
    return(NULL)
  }

  enrich_df$p.adjust <- as.numeric(enrich_df$p.adjust)
  enrich_df <- enrich_df[is.finite(enrich_df$p.adjust) & !is.na(enrich_df$p.adjust), , drop = FALSE]
  if (nrow(enrich_df) == 0) {
    warning("No valid p.adjust values to plot.")
    return(NULL)
  }

  ggplot2::ggplot(enrich_df, ggplot2::aes(x = .data$Community, y = .data$Description,
                                         size = -log10(.data$p.adjust),
                                         color = .data$Community)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = title,
                  x = "Community",
                  y = paste(ontology, "Terms"),
                  size = "-log10(p.adjust)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

