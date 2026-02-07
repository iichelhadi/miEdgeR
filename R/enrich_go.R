#' Perform Enrichment Analysis for Communities
#'
#' Performs enrichment analysis for each community using specified semantic space.
#'
#' @param large_communities A list of communities OR a list of bin results from detect_communities_across_bins()
#'        where each element contains $bin and $communities.
#' @param semantic_space Character, one of: "GO", "KEGG", "DO", "Reactome" (default: "GO").
#' @param ontology Character, GO ontology: "BP", "CC", "MF", "ALL" (default: "BP").
#' @param pvalueCutoff Numeric, p-value cutoff for enrichment (default: 0.05).
#' @param qvalueCutoff Numeric, q-value cutoff for enrichment (default: 0.2).
#' @param verbose Logical, whether to message progress (default: FALSE).
#' @return A named list of enrichment result objects (clusterProfiler/DOSE/ReactomePA).
#' @export
enrich_go <- function(large_communities,
                      semantic_space = "GO",
                      ontology = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      verbose = FALSE) {

  valid_ontologies <- c("BP", "CC", "MF", "ALL")
  if (!ontology %in% valid_ontologies) {
    stop("ontology must be one of: ", paste(valid_ontologies, collapse = ", "))
  }

  valid_spaces <- c("GO", "KEGG", "DO", "Reactome")
  if (!semantic_space %in% valid_spaces) {
    stop("semantic_space must be one of: ", paste(valid_spaces, collapse = ", "))
  }

  if (!is.list(large_communities)) stop("large_communities must be a list.")

  .n_terms <- function(x) {
    if (is.null(x)) return(0L)
    if (!("result" %in% slotNames(x))) return(0L)
    nrow(x@result)
  }

  .filter_genes <- function(genes) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    filter_housekeeping(genes)  # internal package function; do NOT use miEdgeR:: inside package
  }

  .symbol_to_entrez <- function(symbols) {
    # Return unique ENTREZIDs; NULL if none map
    gm <- suppressMessages(
      clusterProfiler::bitr(
        symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    )
    if (is.null(gm) || nrow(gm) == 0) return(NULL)
    unique(gm$ENTREZID)
  }

  .do_enrich <- function(filtered_genes) {
    if (length(filtered_genes) < 5) return(NULL)

    if (semantic_space == "GO") {
      er <- suppressMessages(
        clusterProfiler::enrichGO(
          gene = filtered_genes,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = ontology,
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      )
      if (.n_terms(er) > 1) er <- suppressMessages(clusterProfiler::simplify(er))
      return(er)
    }

    entrez <- .symbol_to_entrez(filtered_genes)
    if (is.null(entrez) || length(entrez) < 5) return(NULL)

    if (semantic_space == "KEGG") {
      return(suppressMessages(
        clusterProfiler::enrichKEGG(
          gene = entrez,
          organism = "hsa",
          keyType = "kegg",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      ))
    }

    if (semantic_space == "DO") {
      return(suppressMessages(
        DOSE::enrichDO(
          gene = entrez,
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      ))
    }

    # Reactome
    suppressMessages(
      ReactomePA::enrichPathway(
        gene = entrez,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    )
  }

  enrich_results <- list()

  is_bin_results <- all(vapply(
    large_communities,
    function(x) is.list(x) && !is.null(x$communities) && !is.null(x$bin),
    logical(1)
  ))

  if (is_bin_results) {
    for (bin_result in large_communities) {
      bin <- as.character(bin_result$bin)
      comms <- bin_result$communities
      if (!is.list(comms) || length(comms) == 0) next

      for (i in seq_along(comms)) {
        genes <- comms[[i]]
        filtered_genes <- .filter_genes(genes)

        er <- .do_enrich(filtered_genes)
        if (.n_terms(er) > 0) {
          enrich_results[[paste(bin, "Module", i, sep = "_")]] <- er
        } else if (verbose) {
          message("No enrichment: bin=", bin, " module=", i, " (genes=", length(filtered_genes), ")")
        }
      }
    }
  } else {
    for (i in seq_along(large_communities)) {
      genes <- large_communities[[i]]
      filtered_genes <- .filter_genes(genes)

      er <- .do_enrich(filtered_genes)
      if (.n_terms(er) > 0) {
        enrich_results[[paste0("Module_", i)]] <- er
      } else if (verbose) {
        message("No enrichment: module=", i, " (genes=", length(filtered_genes), ")")
      }
    }
  }

  enrich_results
}

