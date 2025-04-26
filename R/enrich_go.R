#' Perform Enrichment Analysis for Communities
#'
#' Performs enrichment analysis for each community using specified semantic space.
#'
#' @param large_communities A list of communities or a list of community detection results (from detect_communities_across_bins).
#' @param semantic_space Character, semantic space to use ("GO", "KEGG", "DO", "Reactome") (default: "GO").
#' @param ontology Character, GO ontology to use ("BP", "CC", "MF", or "ALL") (default: "BP").
#' @param pvalueCutoff Numeric, p-value cutoff for enrichment (default: 0.05).
#' @param qvalueCutoff Numeric, q-value cutoff for enrichment (default: 0.2).
#' @return A named list of enrichment results.
#' @export
#' @importFrom clusterProfiler enrichGO simplify enrichKEGG bitr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom DOSE enrichDO
#' @importFrom ReactomePA enrichPathway
enrich_go <- function(large_communities,
                      semantic_space = "GO",
                      ontology = "BP",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  # Validate ontology
  valid_ontologies <- c("BP", "CC", "MF", "ALL")
  if (!ontology %in% valid_ontologies) {
    stop("Ontology must be one of: ", paste(valid_ontologies, collapse = ", "))
  }

  # Validate semantic_space
  valid_spaces <- c("GO", "KEGG", "DO", "Reactome")
  if (!semantic_space %in% valid_spaces) {
    stop("Semantic space must be one of: ", paste(valid_spaces, collapse = ", "))
  }

  enrich_results <- list()

  # Check if large_communities is a list of community results (e.g., from detect_communities_across_bins)
  if (all(sapply(large_communities, function(x) "communities" %in% names(x)))) {
    # Process each bin's communities
    for (bin_result in large_communities) {
      bin <- bin_result$bin
      comms <- bin_result$communities
      for (i in seq_along(comms)) {
        genes <- comms[[i]]
        filtered_genes <- miEdgeR::filter_housekeeping(genes)
        cat("Community", i, "in bin", bin, "- Genes after filtering:", length(filtered_genes), "\n")

        # Perform enrichment based on semantic_space
        if (semantic_space == "GO") {
          enrich_result <- enrichGO(
            gene = filtered_genes,
            OrgDb = org.Hs.eg.db,
            keyType = "SYMBOL",
            ont = ontology,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          )
          if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
            enrich_result <- simplify(enrich_result)
          }
        } else if (semantic_space == "KEGG") {
          # Convert gene symbols to Entrez IDs
          gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          if (nrow(gene_map) == 0) {
            cat("No Entrez IDs mapped for Community", i, "in bin", bin, "\n")
            next
          }
          enrich_result <- enrichKEGG(
            gene = gene_map$ENTREZID,
            organism = "hsa",
            keyType = "ENTREZID",
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          )
        } else if (semantic_space == "DO") {
          # Convert gene symbols to Entrez IDs
          gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          if (nrow(gene_map) == 0) {
            cat("No Entrez IDs mapped for Community", i, "in bin", bin, "\n")
            next
          }
          enrich_result <- DOSE::enrichDO(
            gene = gene_map$ENTREZID,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          )
        } else if (semantic_space == "Reactome") {
          # Convert gene symbols to Entrez IDs
          gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          if (nrow(gene_map) == 0) {
            cat("No Entrez IDs mapped for Community", i, "in bin", bin, "\n")
            next
          }
          enrich_result <- ReactomePA::enrichPathway(
            gene = gene_map$ENTREZID,
            organism = "human",
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          )
        }

        if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
          enrich_results[[paste(bin, "Group", i, sep = "_")]] <- enrich_result
        } else {
          cat("No significant enrichment for Community", i, "in bin", bin, "\n")
        }
      }
    }
  } else {
    # Process as a simple list of communities
    for (i in seq_along(large_communities)) {
      genes <- large_communities[[i]]
      filtered_genes <- miEdgeR::filter_housekeeping(genes)
      cat("Community", i, "- Genes after filtering:", length(filtered_genes), "\n")

      if (semantic_space == "GO") {
        enrich_result <- enrichGO(
          gene = filtered_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = ontology,
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
        if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
          enrich_result <- simplify(enrich_result)
        }
      } else if (semantic_space == "KEGG") {
        # Convert gene symbols to Entrez IDs
        gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(gene_map) == 0) {
          cat("No Entrez IDs mapped for Community", i, "\n")
          next
        }
        enrich_result <- enrichKEGG(
          gene = gene_map$ENTREZID,
          organism = "hsa",
          keyType = "ENTREZID",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      } else if (semantic_space == "DO") {
        # Convert gene symbols to Entrez IDs
        gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(gene_map) == 0) {
          cat("No Entrez IDs mapped for Community", i, "\n")
          next
        }
        enrich_result <- DOSE::enrichDO(
          gene = gene_map$ENTREZID,
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      } else if (semantic_space == "Reactome") {
        # Convert gene symbols to Entrez IDs
        gene_map <- clusterProfiler::bitr(filtered_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(gene_map) == 0) {
          cat("No Entrez IDs mapped for Community", i, "\n")
          next
        }
        enrich_result <- ReactomePA::enrichPathway(
          gene = gene_map$ENTREZID,
          organism = "human",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )
      }

      if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
        enrich_results[[paste("Community", i)]] <- enrich_result
      } else {
        cat("No significant enrichment for Community", i, "\n")
      }
    }
  }

  # Print enrichment results
  for (name in names(enrich_results)) {
    cat("Enrichment for", name, ":\n")
    print(head(enrich_results[[name]]@result, 4))
  }

  enrich_results
}
