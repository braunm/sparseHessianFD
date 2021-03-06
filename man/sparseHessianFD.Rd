% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseHessianFD-class.R
\docType{class}
\name{sparseHessianFD}
\alias{sparseHessianFD}
\title{sparseHessianFD}
\description{
A reference class for computing sparse Hessians
}
\details{
The sparseHessianFD function calls the initializer for the
sparseHessianFD class, and returns a sparseHessianFD object.
\preformatted{sparseHessianFD(x, fn, gr, rows, cols, delta, index1, complex, ...)}
The function, gradient and sparsity pattern are
declared as part of the initialization.

Once initialized, the $hessian method will evaluate the Hessian at x.
\preformatted{
obj <- sparseHessian(x, fn, gr, rows, cols, ...)
obj$hessian(x)
}
For convenience, the class provides wrapper methods to the \code{fn} and \code{gr} functions that were specified in the initializer.
\preformatted{
obj$fn(P) ## wrapper to objective function
obj$gr(P) ## wrapper to gradient
obj$fngr(P) ## list of obj function and gradient
obj$fngrhs(P) ## list of obj function, gradient and Hessian.
}
\subsection{Arguments to initializer}{
\describe{
\item{x}{an vector at which the function, gradient and Hessian are initialized and tested.}
\item{fn, gr}{R functions that return the function value and
gradient, evaluated at x.}
\item{rows, cols}{Numeric vectors with row and column indices of
the non-zero elements in the lower triangle (including diagonal) of
the Hessian.}
\item{delta}{The perturbation amount for finite difference (or
complex step) of the gradient to compute the Hessian.  Defaults to 1e-07.}
\item{index1}{TRUE if rows and cols use 1-based (R format)
indexing (FALSE for 0-based (C format) indexing.}
\item{complex}{TRUE if Hessian will be computed using the complex
step method, and FALSE (default) if using finite differences. If
TRUE, both fn and gr must accept complex arguments and return complex values.}
\item{...}{other arguments to be passed to fn and gr.}
}
}
Other methods are described below.  Do not access any of the
    fields directly.  The internal structure is subject to change in future versions.
}
\section{Fields}{

\describe{
\item{\code{fn1}}{A closure for calling fn(x, ...).}

\item{\code{gr1}}{A closure for calling gr(x, ...).}

\item{\code{iRow,jCol}}{Numeric vectors with row and column indices of
the non-zero elements in the lower triangle (including diagonal) of
the Hessian.}

\item{\code{delta}}{The perturbation amount for finite differencing of
the gradient to compute the Hessian. Defaults to 1e-07.}

\item{\code{index1}}{TRUE if rows and cols use 1-based (R format)
indexing (FALSE for 0-based (C format) indexing.}

\item{\code{complex}}{TRUE if Hessian will be computed using the complex
step method, and FALSE (default) if using finite differences.}

\item{\code{D}}{raw finite differences (internal use only)}

\item{\code{nvars}}{Number of variables (length of x)}

\item{\code{nnz}}{Number of non-zero elements in the lower triangle of the Hessian.}

\item{\code{ready}}{TRUE if object has been initialized, and Hessian has
been partitioned.}

\item{\code{idx,pntr}}{Column indices and row pointers for non-zero
elements in lower triangle of the permuted Hessian.  Row-oriented
compressed storage.}

\item{\code{colors}}{A vector representation of the partitioning of the columns.
There are nvars elements, one for each column of the permuted
Hessian.  The value corresponds to the "color" for that column.}

\item{\code{perm,invperm}}{Permutation vector and its inverse}
}}

\section{Methods}{

\describe{
\item{\code{fn(x)}}{Return function value, evaluated at x: fn(x, ...)}

\item{\code{fngr(x)}}{Return list of function value and gradient, evaluated at x}

\item{\code{fngrhs(x)}}{Return list of function value, gradient, and Hessian, evaluated at x}

\item{\code{get_invperm()}}{Return integer vector of inverse of permutation used for computing Hessian}

\item{\code{get_nnz()}}{Return number of non-zero elements in lower triangle of Hessian}

\item{\code{get_nvars()}}{Return dimension (number of rows or columns) of Hessian}

\item{\code{get_pattern()}}{Return pattern matrix of lower triangle of Hessian}

\item{\code{get_perm()}}{Return integer vector of permutation used for computing Hessian}

\item{\code{get_perm_pattern()}}{Return pattern matrix of lower triangle of *permuted* Hessian}

\item{\code{gr(x)}}{Return gradient, evaluated at x:  gr(x,...)}

\item{\code{hessian(x)}}{Return sparse Hessian, evaluated at x, as a dgCMatrix object.}

\item{\code{initialize(
  x,
  fn,
  gr,
  rows,
  cols,
  delta = 1e-07,
  index1 = TRUE,
  complex = FALSE,
  ...
)}}{Initialize object with functions to compute the objective function and gradient (fn and gr), row and column indices of non-zero elements (rows and cols), an initial variable vector x at which fn and gr can be evaluated, a finite differencing parameter delta, flags for 0 or 1-based indexing (index1), whether the complex step method will be used, and other arguments (...) to be passed to fn and gr.}

\item{\code{partition()}}{Return the partitioning used to compute finite differences}

\item{\code{pointers(out.index1 = index1)}}{Return list with indices (idx) and pointers (pntr) for sparsity pattern of the compressed sparse Hessian.  Since the Hessian is symmetric, the indices and pointers for row-oriented and column-oriented storage patterns are the same.}
}}

\examples{
## Log posterior density of hierarchical binary choice model. See vignette.
set.seed(123)
data("binary_small")
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
T <- binary[["T"]]
P <- rnorm((N+1)*k)
priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))
true.hess <- binary.hess(P, binary, priors)
pattern <- Matrix.to.Coord(Matrix::tril(true.hess))
str(pattern)
obj <- sparseHessianFD(P, fn=binary.f, gr=binary.grad,
       rows=pattern[["rows"]], cols=pattern[["cols"]],
                      data=binary, priors=priors)
hs <- obj$hessian(P)
all.equal(hs, true.hess)

f <- obj$fn(P) ## obj function
df <- obj$gr(P) ## gradient
fdf <- obj$fngr(P) ## list of obj function and gradient
fdfhs <- obj$fngrhs(P) ## list of obj function, gradient and Hessian.
}
