% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mi_network.R
\name{get_expr_data}
\alias{get_expr_data}
\title{Get Expression Data for a Cluster}
\usage{
get_expr_data(
  seurat_obj,
  cluster_id,
  cluster_field,
  assay_name = "RNA",
  counts_layer = "counts",
  data_layer = "data",
  top_n_genes = 5000,
  max_cells = 1500,
  min_expr_pct = 0.05
)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{cluster_id}{Cluster ID}

\item{cluster_field}{Metadata column name for clusters}

\item{assay_name}{Assay name (default: "RNA")}

\item{counts_layer}{Layer for counts data (default: "counts")}

\item{data_layer}{Layer for normalized data (default: "data")}

\item{top_n_genes}{Number of top variable genes}

\item{max_cells}{Maximum number of cells}

\item{min_expr_pct}{Minimum expression percentage}
}
\value{
Expression matrix
}
\description{
Get Expression Data for a Cluster
}
