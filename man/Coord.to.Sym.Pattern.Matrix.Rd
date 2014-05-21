\name{Coord.to.Sym.Pattern.Matrix}
\alias{Coord.to.Sym.Pattern.Matrix}
\title{Converts row and column indices to a symmetric pattern Matrix object of Matrix class}
\description{
Converts row and column indices to a symmetric pattern Matrix object of Matrix class
	}
\usage{
Coord.to.Sym.Pattern.Matrix(H, nvars)
}

\arguments{
	\item{H}{a list containing matrix structure data.  See details.}
	\item{nvars}{the number of rows (and columns) in the matrix.}
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
Only the non-zero elements in the lower triangle should be included.
}




 


