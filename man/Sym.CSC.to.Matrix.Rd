\name{Sym.CSC.to.Matrix}
\alias{Sym.CSC.to.Matrix}
\title{Converts symmetric Hessian stored in CSC format to an object of Matrix class}
\description{
Converts symmetric Hessian stored in CSC format to an object of Matrix class.
	}
\usage{
Sym.CSC.to.Matrix(H, nvars)
}

\arguments{
	\item{H}{a list containing Hessian data.  See details.}
	\item{nvars}{the number of rows (and columns) in the matrix.
    }
}
\value{
An object of Matrix class.
}
\details{
H is a list with three elements.
\describe{
\item{iRow}{Row indices for each of the non-zero elements in the lower triangle of H}
\item{jpntr}{A vector of length \code{nvars+1}.  jpntr[j] is the index of the element in vals that is the first non-zero element in the jth column of the matrix.}
\item{vals}{The values of the non-zero elements in the lower triangle of the matrix.}
}
}




 


