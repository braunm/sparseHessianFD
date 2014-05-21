\name{Matrix.to.Coord}
\alias{Matrix.to.Coord}
\title{Takes an object of class Matrix and returns the row and column indices of the non-zero elements in the lower triangle.}
\description{
Takes an object of class Matrix and returns the row and column indices of the non-zero elements in the lower triangle
	}
\usage{
Matrix.to.Coord(M)
}

\arguments{
	\item{M}{a symmetric Matrix object.}
    }
\value{
A list suitable for using as the \code{hess.struct} argument in \code{get.sparse.hessian}.  There are three elements:
\item{iRow}{Row indices for each of the non-zero elements in the lower triangle of \code{M}. }
\item{jpntr}{A vector of length \code{nvars+1}.  jpntr[j] is the index of the element in vals that is the first non-zero element in the jth column of \code{M}.}
\item{vals}{The values of the non-zero elements in the lower triangle of \code{M}.}
}




 


