% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{RfastTOMdist}
\alias{RfastTOMdist}
\title{RfastTOMdist}
\usage{
RfastTOMdist(A, mat_mult_method = c("Rfast", "RcppEigen"))
}
\arguments{
\item{A}{A N x N adjacency matrix where N is the number of genes. values range from -1:1 with higher values indicating highly similar genes. Often correlation ^exponent, but could be angular distance, mutual information or Euclidean distance}

\item{mat_mult_method}{method for large matrix multiplication, "Rfast" (default) or "RcppEigen" (see \code{details} in \code{\link[=icwgcna]{icwgcna()}})}
}
\value{
A N x N distance matrix with smaller values indicating more related genes.
}
\description{
distance based on the topological overlap map from \href{https://pubmed.ncbi.nlm.nih.gov/12202830/}{Ravasz, E., Somera, A., Mongru, D., Oltvai, Z. and Barab´asi, A. (2002). Science}
Implemented to using the \link[Rfast:Rfast-package]{Rfast} functions to  speed things up since we will be computing this up to 25 times
}
