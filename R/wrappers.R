## wrappers.R --   Part of the sparsehessianFD package for the R programming language.
##
## Copyright (C) 2015 Michael Braun


#' @name sparseHessianFD.new
#' @title Create and initialize a new sparseHessianFD object
#' @param x vector of variables at which to evaluate value, gradient and Hessian during initialization.
#' @param fn R function that returns function value
#' @param gr R function that returns the gradient of the function
#' @param rows Integer vector of row indices of non-zero elements of the lower triangle of the Hessian
#' @param cols Integer vector of column indices of non-zero elements of the lower triangle of the Hessian
#' @param direct If TRUE, use direct method.  Otherwise, use indirect/substitution method.
#' @param eps Difference parameter for finite differencing
#' @param ... Other parameters to be passed to fn and gr.
#' @return An object of class sparseHessianFD
#' @details The function new.sparse.hessian.obj has been deprecated.  Use sparseHessianFD.new instead.
#' @seealso sparseHessianFD-class
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
