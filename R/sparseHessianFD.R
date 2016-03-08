#' @name sparseHessianFD-package
#' @aliases sparseHessianFD-package
#' @docType package
#' @title Estimate sparse Hessians using finite differences of
#' gradients.
#' @description Estimate sparse Hessian
#' @details This package contains methods to return the Hessian of
#' a function in sparse Matrix format.  The user must supply the
#' objective function, the gradient, and the row and column indices of
#' the non-zero elements of the lower triangle of the Hessian (i.e.,
#' the sparsity structure must be known in advance).

#' @references
#' Coleman, Thomas F, and Jin-Yi Cai. 1986.  The Cyclic Coloring Problem
#' and Estimation of Sparse Hessian Matrices.  SIAM Journal on Algebraic
#' Discrete Methods 7 (2): 221-235
#'
#' Coleman, Thomas F, Burton S Garbow, and Jorge J Mor\'e. 1985. Software
#' for Estimating Sparse Hessian Matrices. ACM Transaction on
#' Mathematical Software 11 (4) (December): 363-377.
#'
#' Coleman, Thomas F, Burton S Garbow, and Jorge J Mor\'e. 1985. Algorithm
#' 636:  FORTRAN Subroutines for Estimating Sparse Hessian Matrices. ACM
#' Transactions on Mathematical Software 11 (4): 378.
#'
#' Coleman, Thomas F and Jorge J More\'e.  Estimation of Sparse Hessian
#' Matrices and Graph Coloring Problems.  Mathematical Programming
#' 28 (3) (October): 243-270
#'
#' Powell, M J D and Ph L Toint. 1979. On the Estimation of Sparse
#' Hessian Matrices.  SIAM Journal on Numerical Analysis 16 (6)
#' (December): 1060-1074.
#'
#' @keywords package
#'
#' @useDynLib sparseHessianFD
#' @encoding UTF-8
#' @import Rcpp
#' @import methods
#' @import Matrix
NULL
