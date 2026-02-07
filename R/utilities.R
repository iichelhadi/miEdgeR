#' Filter Housekeeping Genes
#'
#' Removes common housekeeping / technical genes (ribosomal, GAPDH/ACTB, mitochondrial).
#'
#' @param gene_list Character vector of gene names.
#' @return Character vector of filtered gene names.
#' @export
filter_housekeeping <- function(gene_list) {
  gene_list <- as.character(gene_list)
  gene_list[!grepl("^(RPL|RPS|GAPDH|ACTB|MT-)", gene_list)]
}

#' Choose Adaptive Number of Bins
#'
#' Heuristic bin choice (Freedman–Diaconis-like), constrained by min/max bins and
#' a minimum expected cells-per-bin.
#'
#' @param expr_data Numeric vector OR matrix (genes x cells) used to estimate bins.
#' @param min_bins Minimum number of bins.
#' @param max_bins Maximum number of bins.
#' @param min_cells_per_bin Minimum cells per bin.
#' @return Integer number of bins.
#' @export
#' @importFrom stats sd
choose_adaptive_nbins <- function(expr_data, min_bins = 5, max_bins = 30, min_cells_per_bin = 75) {
  if (is.null(dim(expr_data))) n_total <- length(expr_data) else n_total <- ncol(expr_data)
  non_na_vals <- expr_data[!is.na(expr_data)]
  n <- length(non_na_vals)

  if (n < 2 || all(diff(sort(non_na_vals)) == 0)) return(min_bins)

  sd_val <- stats::sd(expr_data, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) return(min_bins)

  h <- 3.5 * sd_val / (n_total^(1/3))
  range_data <- diff(range(expr_data, na.rm = TRUE))
  if (is.na(h) || h == 0) return(min_bins)

  nbins <- ceiling(range_data / h)
  if (is.na(nbins) || nbins <= 0 || is.infinite(nbins)) return(min_bins)

  avg_cells_per_bin <- n_total / nbins
  while (avg_cells_per_bin < min_cells_per_bin && nbins > min_bins) {
    nbins <- nbins - 1
    avg_cells_per_bin <- n_total / nbins
  }

  max(min_bins, min(max_bins, nbins))
}

#' Ensure Graph Has Finite Edge Weights
#'
#' Adds/repairs `E(g)$weight` using `mean_weight` if present, else defaults to 1.
#'
#' @param g An igraph object.
#' @return igraph object with finite numeric `weight` edge attribute.
#' @export
#' @importFrom igraph edge_attr_names E
ensure_weight <- function(g) {
  if (!inherits(g, "igraph")) stop("g must be an igraph object.")
  if (!("weight" %in% igraph::edge_attr_names(g))) {
    if ("mean_weight" %in% igraph::edge_attr_names(g)) {
      igraph::E(g)$weight <- igraph::E(g)$mean_weight
    } else {
      igraph::E(g)$weight <- 1
    }
  }
  w <- igraph::E(g)$weight
  bad <- !is.finite(w) | is.na(w)
  if (any(bad)) igraph::E(g)$weight[bad] <- 1
  g
}

#' Degree as Numeric Vector
#'
#' @param g An igraph object.
#' @return Numeric vector of degrees.
#' @export
#' @importFrom igraph degree
degree_num <- function(g) {
  as.numeric(unlist(igraph::degree(g)))
}

#' Get Assay Matrix (Seurat v4/v5 Compatible)
#'
#' Uses `layer=` when available (SeuratObject v5), else falls back to `slot=`.
#'
#' @param seurat_obj Seurat object.
#' @param assay Assay name.
#' @param layer_or_slot Layer (v5) or slot (v4): e.g. "counts", "data", "scale.data".
#' @return Matrix-like object.
#' @export
#' @importFrom Seurat GetAssayData
get_assay_mat <- function(seurat_obj, assay, layer_or_slot) {
  fn <- Seurat::GetAssayData
  if ("layer" %in% names(formals(fn))) {
    return(Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot))
  }
  Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
}

#' Canonicalize Undirected Edge Data Frame
#'
#' Forces undirected edge representation by sorting endpoints within each row.
#'
#' @param df Data frame with columns `from` and `to` (and optionally `weight`).
#' @return Data frame with canonicalized `from`/`to`.
#' @export
canonical_edge_df <- function(df) {
  if (!all(c("from", "to") %in% colnames(df))) {
    stop("df must contain columns: from, to")
  }
  out <- df
  out$from <- pmin(out$from, out$to)
  out$to   <- pmax(df$from, df$to)
  out
}

