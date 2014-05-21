\encoding{utf-8}
\name{sparseHessianFD-package}
\alias{sparseHessianFD-package}
\alias{sparseHessianFD}
\docType{package}
\title{
Given function, gradient, and sparsity structure, returns Hessian in sparse matrix format.
}
\description{
This package contains methods to return the Hessian of a function
in sparse Matrix format.  The user must supply the objective function, the 
gradient, and the row and column indices of the non-zero elements of the lower
triangle of the Hessian (i.e., the sparsity structure must be known in advance).
Sparse Hessians occur in many applications, such as log posterior densities of 
hierarchical models under conditional independence assumptions.  This package is
intended to be useful when optimizing objective functions with this structure, using
optimizers than can exploit this sparsity, such as the trustOptim package.

License details are available in the LICENSE file in the source code.

}
\details{
\tabular{ll}{
Package: \tab sparseHessianFD\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2012-11-07\cr
License: \tab See the LICENSE file included in the source code.\cr
LazyLoad: \tab yes\cr
}

}
\author{
Michael Braun <braunm@mit.edu>
}

\references{
Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Software
for Estimating Sparse Hessian Matrices. ACM Transaction on
Mathematical Software 11 (4) (December): 363–377. 

Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Algorithm
636:  FORTRAN Subroutines for Estimating Sparse Hessian Matrices. ACM
Transactions on Mathematical Software 11 (4): 378.
}

\keyword{ package }

