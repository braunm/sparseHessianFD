\name{Coord.to.Pattern.Matrix}
\alias{Coord.to.Pattern.Matrix}
\title{Converts row and column indices to a pattern Matrix object of Matrix class}
\description{
Converts row and column indices to a pattern Matrix object of Matrix class
	}
\usage{
Coord.to.Pattern.Matrix(H, nrows, ncols=nrows)
}

\arguments{
	\item{H}{a list containing matrix structure data.  See details.}
	\item{nrows}{the number of rows in the matrix.}
	\item{ncols}{the number of columns in the matrix.}
    }

\value{
An object of Matrix class.
}
\details{
H is a list with two vectors.
\describe{
\item{iRow}{Row indices for each of the non-zero elements in the matrix.}
\item{jCol}{Column indices for each of the non-zero elements in the matrix.}
}
}




 


