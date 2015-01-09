## wrappers.R --   Part of the sparsehessianFD package for the R programming language.
##
## Copyright (C) 2015 Michael Braun


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
#' @details  Indexing starts at 1, and should include the diagonal
#' elements, even though it is obvious that the diagonal elements are
#' non-zero.  Entries must be ordered by column, and by row within
#' columns.  Do not include any entries for the upper triangle.
#' The algorithms used for estimating sparse Hessians using finite
#' differencing are described in the reference below.
#'
#' This method involves a partitioning and permutation of the Hessian
#' matrix to reduce the number of distinct finite differencing
#' operations of the gradient.  There are two methods for computing
#' the partition and permutation:  direct and indirect.  The direct
#' method tends to require more computation, but may be more
#' accurate.  We recommend the indirect method to start, so the
#' default value is 0.
#'
#' The function new.sparse.hessian.obj has been deprecated.  Use sparseHessianFD.new instead.
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
    if (any(cols > rows)) stop("sparseHessianFD: Some elements in upper triangle.  Provide lower triangle only.")

    obj$hessian.init(rows, cols, as.integer(direct), eps)
    
    return(obj)
    
}


#' @rdname sparseHessianFD.new
new.sparse.hessian.obj <- function(x, fn, gr, hs, fd.method=0L, eps=sqrt(.Machine$double.eps),...) {
    .Deprecated("sparseHessianFD.new")
    if (is.null(hs))
        stop("sparseHessianFD: you must supply structure of the Hessian.")
    if (!is.list(hs))
        stop ("sparseHessianFD:  hs must be a list")
    if (!all.equal(names(hs),c("iRow","jCol")))
        stop ("sparseHessianFD:  Names of hs must be iRow and jCol")
    direct <- as.logical(fd.method)        
    return(sparseHessianFD.new(x, fn, gr, hs$iRows, hs$jCols, direct, eps, ...))
}
