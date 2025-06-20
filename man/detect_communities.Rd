% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mi_network.R
\name{detect_communities}
\alias{detect_communities}
\title{Detect Communities
Detects communities using label propagation and neighborhood overlap.}
\usage{
detect_communities(graph, min_size = 10, overlap_threshold = 0.1)
}
\arguments{
\item{graph}{Input igraph object.}

\item{min_size}{Minimum community size.}

\item{overlap_threshold}{Overlap threshold for adding genes (default: 0.1).}
}
\value{
List with communities and their sizes.
}
\description{
Detect Communities
Detects communities using label propagation and neighborhood overlap.
}
