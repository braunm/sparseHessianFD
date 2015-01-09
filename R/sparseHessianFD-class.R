## sparseHessianFD-class.R -- Part of the sparseHessianFD package 
## Copyright (C) 2014-2015 Michael Braun
## See LICENSE file for details.




##' @title sparseHessianFD
##' @name sparseHessianFD-class
##' @description This class provides methods to estimate the Hessian of an objective function.
##' @aliases Rcpp_sparseHessianFD-class
##' @docType class
##' @details
##' Methods to estimate the Hessian of an
##' objective function using finite ifferencing of gradients. An
##' object is created with the names of R functions that return the
##' value and the gradient, and initialized with the row and column
##' indices of the non-zero elements in the lower triangle of the
##' Hessian.
##' The class contains the following methods.  See vignette for usage.
##' \describe{
##' \item{$fn(x)}{\code{signature(x="numeric")}: returns fn(x)}
##' \item{$gr(x)}{\code{signature(x="numeric")}: returns gr(x)}
##' \item{$fngr(x)}{\code{signature(x="numeric")}: returns list of fn(x) and gr(x)}
##' \item{$hessian(x)}{\code{signature(x="numeric")}: returns sparse
##' Hessian as dgCMatrix object}
##' \item{$hessian.init(rows, cols, direct, eps)}{ Used internally to initialize Hessian with sparsity pattern}
##' \item{$nnz()}{ Number of non-zero elements in lower triangle of Hessian, as provided by sparsity pattern}
##' \item{$nvars()}{ Length of parameter vector that was provided to constructor}
##' }
##' @export sparseHessianFD
sparseHessianFD <- setRcppClass("sparseHessianFD", module="sparseHessianFD")
