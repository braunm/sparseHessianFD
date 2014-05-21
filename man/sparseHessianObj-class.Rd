\name{sparseHessianObj-class}
\alias{sparseHessianObj-class}
\alias{$,sparseHessianObj-method}
\alias{new, sparseHessianObj-method}
\title{Class to store function, gradient and Hessian structure information}
\docType{class}
\description{
	The class \code{"sparseHessianObj"} contains an external pointer to a C++ class that contains the methods that are used to compute the sparse Hessian.  It can be constructed directly using \code{"new"}, or through the \code{"new.sparse.hessian.obj"} function.  The latter approach includes some additional argument checking. 
	}

\section{Slots}{
	\describe{
	\item{\code{pointer}:}{Object of class \code{"Rcpp:::externalptr"}}
	}
}

\section{Constructing objects}{
\describe{
Objects of the class can be constructed directly by:

\item{}{\code{new("sparseHessianObj", nvars, fn, gr, hs.struct, fd.method, eps)}}

\item{nvars}{The length of the function argument vector x}
\item{fn}{A function that takes a single argument and returns the value of the objective function, evaluated at x.}
\item{gr}{A function that takes a single argument and returns the gradient of the objective function, evaluated at x.}
\item{hs.struct}{A list of two integer vectors that describe the sparsity structure of the Hessian. See details.}
\item{fd.method}{An integer, 0 for the indirect method and 1 for the direct method, as described in Coleman, et.al (1985). The indirect method is somewhat faster, and suitable for most applications.}
\item{eps}{The perturbation about for finite differencing.  Typically \code{sqrt(.Machine$double.eps)}.}
}
}

\seealso{new.sparse.hessian.obj}



\section{Methods}{
\describe{
\item{$fn(x)}{\code{signature(x="sparseHessianObj")}: returns fn(x)}
\item{$gr(x)}{\code{signature(x="sparseHessianObj")}: returns gr(x)}
\item{$fngr(x)}{\code{signature(x="sparseHessianObj")}: returns list of fn(x) and gr(x)}
\item{$hessian(x)}{\code{signature(x="sparseHessianObj")}: returns hessian(x) as \code{"dgCMatrix"}}
\item{$all(x)}{\code{signature(x="sparseHessianObj")}: returns list of fn(x), gr(x) and hessian}
}
}




