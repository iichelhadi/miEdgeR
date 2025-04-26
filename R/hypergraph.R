#' Create Hypergraph from Gene Groups
#' @param group_list List of gene groups
#' @param all_genes All genes
#' @return List containing hypergraph, incidence matrix, genes, and groups
#' @export
make_hypergraph <- function(group_list, all_genes) {
  all_genes <- unique(all_genes)
  hyperedges <- lapply(group_list, function(group) {
    valid_genes <- intersect(group, all_genes)
    if (length(valid_genes) == 0) return(NULL)
    hypergraph::Hyperedge(valid_genes)
  })
  hyperedges <- Filter(Negate(is.null), hyperedges)
  hg <- if (length(hyperedges) > 0) hypergraph::Hypergraph(all_genes, hyperedges) else hypergraph::Hypergraph(all_genes, list())
  
  incidence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(group_list))
  rownames(incidence_matrix) <- all_genes
  colnames(incidence_matrix) <- names(group_list)
  for (i in seq_along(group_list)) {
    valid_genes <- intersect(group_list[[i]], all_genes)
    if (length(valid_genes) > 0) {
      incidence_matrix[valid_genes, i] <- 1
    }
  }
  
  list(
    hypergraph = hg,
    incidence_matrix = incidence_matrix,
    genes = all_genes,
    groups = group_list
  )
}
