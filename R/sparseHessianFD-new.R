## matrices.R -- Part of the sparseHessianFD package 
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.


#' @name sparseHessianFD.new
#' @title Create and initialize a new sparseHessianFD object
#' @param x A intital vector of variables at which to evaluate value, gradient
#' and Hessian during initialization.
#' @param fn R function that returns function value
#' @param gr R function that returns the gradient of the function
#' @param rows Integer vector of row indices of non-zero elements of
#' the lower triangle of the Hessian
#' @param cols Integer vector of column indices of non-zero elements
#' of the lower triangle of the Hessian 
#' @param direct If TRUE, use direct method for computatation.  Otherwise, use
#' indirect/substitution method.  See references.
#' @param eps The perturbation amount for finite differencing of the gradient to compute the Hessian. Defaults to sqrt(.Machine$double.eps).
#' @param ... Other parameters to be passed to fn and gr.
#' @return An object of class sparseHessianFD
#' @details  Indexing starts at 1, and must include the diagonal
#' elements.  Any upper triangle coordinates will be removed, and
#' replaced with their lower triangle counterparts. 
#' The algorithms used for estimating sparse Hessians using finite
#' differencing are described in the reference below.
#'
#' This method involves a partitioning and permutation of the Hessian
#' matrix to reduce the number of distinct finite differencing
#' operations of the gradient.  There are two methods for computing
#' the partition and permutation:  direct and indirect.  The direct
#' method tends to require more computation, but may be more
#' accurate.  We recommend the indirect method to start, so the
#' default value is FALSE.
#'
#' Here is the description of the two methods, as included in the
#' original ACM TOMS Fortran code:
#'
#'   The direct method (method = 1) first determines a partition of
#'   the columns of symmetric matrix A such that two columns in a
#'   group have a non-zero element in row K only if column K is in an
#'   earlier group. Using this partition, the subroutine then computes
#'   a symmetric permutation of A consistent with the determination of
#'   A by a lower triangular substitution method. 
#' 
#'   The indirect method first computes a symmetric permutation of A
#'   which minimizes the maximum number of non-zero elements in any 
#'   row of L, where L is the lower triangular part of the permuted
#'   matrix. The subroutine then partitions the columns of L into
#'   groups such that columns of L in a group do not have a non-zero
#'   in the same row position.
#'
#'
#' The function new.sparse.hessian.obj has been deprecated.  Use
#' sparseHessianFD.new instead. 
#' @seealso sparseHessianFD-class
#' @references Coleman, Thomas F, Burton S Garbow, and Jorge J More
#' 1985. Software for Estimating Sparse Hessian Matrices. ACM
#' Transaction on Mathematical Software 11 (4) (December): 363-377.
#' @export
sparseHessianFD.new <- function(x, fn, gr, rows, cols, direct=FALSE,
                            eps=sqrt(.Machine$double.eps), ...) {
    
    
    stopifnot(is.function(fn))
    stopifnot(is.function(gr))
    
    fn1 <- function(x) fn(x,...)  ## create closure
    gr1 <- if (!is.null(gr)) function(x) gr(x,...)
    
    ## test fn and gr

    stopifnot(is.finite(fn1(x)))
    gradient <- gr1(x)
    stopifnot(all(is.finite(gradient)))
    stopifnot(length(gradient)==length(x))
    stopifnot(!is.null(rows))
    stopifnot(!is.null(cols))
    
    obj <- new("sparseHessianFD", length(x), fn1, gr1)
    rows <- as.integer(rows)
    cols <- as.integer(cols)
    stopifnot(length(rows)==length(cols))
    stopifnot(all(is.finite(rows)) & all(is.finite(cols)))
    if (any(cols > rows)) {
        warning("sparseHessianFD:  Some elements are in upper triangle and will be deleted.")
        cat("Provide lower triangle only.\n")
        ww <- which(cols > rows)
        rows <- rows[-ww]
        cols <- cols[-ww]
    }
    obj$hessian.init(rows, cols, as.integer(direct), eps)
    
    return(obj)
    
}

