#' @name sparseHessianFD-package
#' @aliases sparseHessianFD-package
#' @docType package
#' @title Estimate sparse Hessians using finite differences of
#' gradients.
#' @details This package contains methods to return the Hessian of
#' a function in sparse Matrix format.  The user must supply the
#' objective function, the gradient, and the row and column indices of
#' the non-zero elements of the lower triangle of the Hessian (i.e.,
#' the sparsity structure must be known in advance). Sparse Hessians
#' occur in many applications, such as log posterior densities of
#' hierarchical models under conditional independence
#' assumptions. This package is intended to be useful when optimizing
#' objective functions with this structure, using optimizers than can
#' exploit this sparsity, such as the trustOptim package.
#'
#' License details are available in the LICENSE file in the source code.
#' @references
#' Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Software
#' for Estimating Sparse Hessian Matrices. ACM Transaction on
#' Mathematical Software 11 (4) (December): 363-377. 
#'
#' Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Algorithm
#' 636:  FORTRAN Subroutines for Estimating Sparse Hessian Matrices. ACM
#' Transactions on Mathematical Software 11 (4): 378.
#' @keywords package
#'
#' @useDynLib sparseHessianFD
#' @encoding UTF-8
#' @import Rcpp
#' @import methods
#' @import Matrix
NULL
