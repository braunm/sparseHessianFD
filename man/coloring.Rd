% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coloring.R
\name{coloring}
\alias{coloring}
\title{Triangular partitioning of variables}
\usage{
coloring(L)
}
\arguments{
\item{L}{sparsity pattern of the Hessian as a lower triangular pattern
matrix}
}
\value{
Integer vector of length nvars with color assignments for each
variable.
}
\description{
cyclic coloring from a lower triangular pattern matrix
}
\details{
For internal use.  Exported in order for replication files for
JSS article to run.
}
