#' Create Hypergraph from Gene Groups
#'
#' @param group_list List of gene groups (character vectors).
#' @param all_genes Character vector of all genes to include as vertices.
#' @return List containing hypergraph, incidence matrix, genes, and groups
#' @export
make_hypergraph <- function(group_list, all_genes) {
  if (!is.list(group_list)) stop("group_list must be a list.")
  all_genes <- unique(as.character(all_genes))
  all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
  if (length(all_genes) == 0) stop("all_genes is empty after cleaning.")

  # Clean + intersect each group, drop empties
  groups_clean <- lapply(group_list, function(g) {
    g <- unique(as.character(g))
    g <- g[!is.na(g) & nzchar(g)]
    intersect(g, all_genes)
  })
  keep <- vapply(groups_clean, function(g) length(g) > 0, logical(1))
  groups_clean <- groups_clean[keep]

  if (length(groups_clean) == 0) {
    hg <- hypergraph::Hypergraph(all_genes, list())
    incidence_matrix <- matrix(0, nrow = length(all_genes), ncol = 0)
    rownames(incidence_matrix) <- all_genes
    return(list(
      hypergraph = hg,
      incidence_matrix = incidence_matrix,
      genes = all_genes,
      groups = list()
    ))
  }

  # Ensure stable names
  gnames <- names(groups_clean)
  if (is.null(gnames) || any(!nzchar(gnames))) {
    names(groups_clean) <- paste0("Module ", seq_along(groups_clean))
  }

  # Build hyperedges and hypergraph (aligned to groups_clean)
  hyperedges <- lapply(groups_clean, function(g) hypergraph::Hyperedge(g))
  hg <- hypergraph::Hypergraph(all_genes, hyperedges)

  # Incidence matrix aligned to groups_clean
  incidence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(groups_clean))
  rownames(incidence_matrix) <- all_genes
  colnames(incidence_matrix) <- names(groups_clean)

  for (i in seq_along(groups_clean)) {
    incidence_matrix[groups_clean[[i]], i] <- 1
  }

  list(
    hypergraph = hg,
    incidence_matrix = incidence_matrix,
    genes = all_genes,
    groups = groups_clean
  )
}

