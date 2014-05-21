\name{new.sparse.hessian.obj}
\alias{new.sparse.hessian.obj}
\title{Construct new object of class sparseHessianObj}
\description{Constructs a new object of class sparseHessianObj that can be used to compute the Hessian of a function, given the gradient function and sparsity structure.}
\usage{
new.sparse.hessian.obj(x, fn, gr, hs, fd.method=0L, eps=sqrt(.Machine$double.eps), ...)
}
\arguments{
	\item{x}{An initial value of the argument to the objective function.  Numeric vector.}
	\item{fn}{A scalar-value function that takes x as its first argument.}
	\item{gr}{A function that returns the gradient of \code{fn(x,...)}, where x is the first argument.}
	\item{hs}{A list of two integer vectors that represents the sparsity structure of the Hessian of fn(x,...).  See details.}
	\item{fd.method}{This flag is 1 if using the direct method for computing a finite-differenced Hessian.  Otherwise, the indirect method is used. Default value is 0. See details.}
	\item{eps}{The perturbation amount for finite differencing of the gradient to compute the Hessian. Defaults to sqrt(.Machine$double.eps).}
	\item{...}{Additional arguments passed to \code{fn} and \code{gr}.}
}

\section{Hessian structure information}{
You must provide information about the sparsity structure of the Hessian. There are two elements to the hs.struct list:
\describe{
	\item{iRow}{row indices of non-zero entries in the \emph{lower triangle} of the Hessian}
	\item{jCol}{column indices of non-zero entries in the \emph{lower triangle} of the Hessian}
}
Indexing starts at 1, and should include the diagonal elements, even though it is obvious that the diagonal elements are non-zero.  Entries must be ordered by column, and by row within columns.  Do not include any entries for the upper triangle.

The algorithms used for estimating sparse Hessians using finite differencing are described in 

Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Software
for Estimating Sparse Hessian Matrices. ACM Transaction on
Mathematical Software 11 (4) (December): 363-377.

This method involves a partitioning and permutation of the Hessian matrix to reduce the number of distinct finite differencing operations of the gradient.  There are two methods for computing the partition and permutation:  direct and indirect.  The direct method tends to require more computation, but may be more accurate.  We recommend the indirect method to start, so the default value is 0.
}
\value{An object of class \code{sparseHessianObj}}
