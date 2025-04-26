#' Filter Housekeeping Genes
#' @param gene_list Vector of gene names
#' @return Filtered gene list excluding housekeeping genes
#' @export
filter_housekeeping <- function(gene_list) {
  gene_list[!grepl("^(RPL|RPS|GAPDH|ACTB)", gene_list)]
}

#' Choose Adaptive Number of Bins
#' @param expr_data Expression data
#' @param min_bins Minimum number of bins
#' @param max_bins Maximum number of bins
#' @param min_cells_per_bin Minimum cells per bin
#' @return Number of bins
#' @export
choose_adaptive_nbins <- function(expr_data, min_bins = 5, max_bins = 30, min_cells_per_bin = 75) {
  if (is.null(dim(expr_data))) n_total <- length(expr_data) else n_total <- ncol(expr_data)
  non_na_vals <- expr_data[!is.na(expr_data)]
  n <- length(non_na_vals)
  if (n < 2 || all(diff(sort(non_na_vals)) == 0)) return(min_bins)
  sd_val <- sd(expr_data, na.rm = TRUE)
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
