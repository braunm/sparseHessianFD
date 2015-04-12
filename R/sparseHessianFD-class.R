## This file is part of the sparseHessianFD package
## Copyright (C) 2015 Michael Braun
##
##' @title sparseHessianFD
##' @name sparseHessianFD-class
##' @description A reference class for computing sparse Hessians
##' @docType class
##' @field fn1 A closure for calling fn(x, ...).
##' @field gr1 A closure for calling gr(x, ...).
##' @field rows,cols Numeric vectors with row and column indices of sparse Hessian. May be just the lower triangle elements or full Hessian, depending on indexLT.  Diagonal elements must be included.
##' @field eps The perturbation amount for finite differencing of the gradient to compute the Hessian. Defaults to sqrt(.Machine$double.eps).
##' @field index1 TRUE if rows and cols use 1-based (R format) indexing (FALSE for 0-based (C format) indexing.  Internal storage is always 0-based
##' @field indexLT TRUE (default) if rows and cols index only the lower triangle of the Hessian.  FALSE if the entire Hessian is indexed.
##' @field D raw finite differences (internal use only)
##' @field nvars Number of variables (length of x)
##' @field nnz Number of non-zero elements in the full Hessian (not just the lower triangle)
##' @field ready TRUE if object has been initialized, and Hessian has been partitioned according to the sparsity pattern.
##' @export sparseHessianFD
sparseHessianFD <-
    setRefClass("sparseHessianFD",
                fields=list(
                    fn1 = "function",
                    gr1 = "function",
                    rows = "numeric",
                    cols = "numeric",
                    eps = "numeric",
                    index1 = "logical",
                    indexLT = "logical",
                    D = "matrix",
                    nvars = "integer",
                    nnz = "integer",
                    ready = "logical",
                    idx = "integer",
                    pntr = "integer",
                    colors = "list",
                    colors_vec = "integer"),
                methods = list(
                    initialize = function(x, fn, gr, rows, cols, direct=NULL,
                                          eps=sqrt(.Machine$double.eps),
                                          index1 = TRUE, indexLT = TRUE, ...) {
                        "Initialize object with functions to compute the objective function and gradient (fn and gr), row and column indices of non-zero elements (rows and cols), an initial variable vector x at which fn and gr can be evaluated, a finite differencing parameter eps, flags for 0 or 1-based indexing (index1), whether sparsity pattern is just for the lower triangle (indexLT), and other arguments (...) to be passed to fn and gr."

                        if (!is.null(direct)) {
                            warning(" 'direct' argument is ignored. Only indirect method is, and will be, supported.")
                        }
                        validate(fn, gr, rows, cols, x, eps, index1, indexLT)

                        initFields(fn1 = function(y) fn(y, ...),
                                   gr1 = function(y) gr(y, ...),
                                   rows = as.integer(rows),
                                   cols = as.integer(cols),
                                   eps = eps,
                                   index1 = index1,
                                   indexLT = indexLT,
                                   nvars = length(x),
                                   nnz = length(rows),
                                   ready = FALSE)


                        ptr <- Coord.to.Pointers(rows, cols, c(nvars, nvars),
                                                 triangle=indexLT,
                                                 lower=TRUE,
                                                 order="symmetric",
                                                 index1=TRUE)

                        idx <<- ptr$idx
                        pntr <<- ptr$pntr

                        colors <<- color_graph(pntr, idx, nvars)
                        color_vec <- rep(0, nvars)
                        for (i in 1:length(colors)) {
                            color_vec[colors[[i]]] <- as.integer(i)
                        }
                        .self$usingMethods("coord2vec")
                        D <<- sapply(colors, coord2vec)
                        ready <<- TRUE
                    },

                    partition = function() {
                        "Return the partitioning used to compute finite differences"
                        if (ready) {
                            res <- colors
                            if (index1) {
                                for (i in 1:length(res)) {
                                    res[[i]] <- colors[[i]]+1
                                }
                            } else {
                                res <- colors
                            }
                        } else {
                            stop("Invalid/incomplete partition")
                            return(NULL)
                        }
                    },

                    validate = function(fn, gr, rows, cols, x,
                                        eps, index1, indexLT, ...) {

                        stopifnot(is.numeric(x),
                                  is.function(fn),
                                  is.function(gr),
                                  is.finite(fn(x, ...)),
                                  !is.null(rows),
                                  !is.null(cols),
                                  length(rows)==length(cols),
                                  all(is.finite(rows)),
                                  all(is.finite(cols)))

                        gradient <- gr(x, ...)

                        I1 <- as.integer(index1)
                        check.index1.row <- min(rows) >= I1 &
                          (max(rows) <= (nvars + I1 - 1))
                        check.index1.col <- min(cols) >= I1 &
                          (max(cols) <= (nvars + I1 - 1))

                        stopifnot(all(is.finite(gradient)),
                                  length(gradient)==nvars,
                                  check.index1.row,
                                  check.index1.col
                                  )

                        if (indexLT & any(cols > rows)) {
                                stop("indexLT == TRUE, but some elements are in upper triangle.")
                            }
                    },

                    fd = function(d, x, grad.x) {
                        gr1(x+d) - grad.x
                    },

                    hessian = function(x) {
                        "Return sparse Hessian, evaluated at x, as a dgCMatrix object."
                        if (ready) {
                            grad.x <- gr1(x)
                            Y <- apply(D, 2, fd, x = x, grad.x = grad.x)
                            res <- subst(Y, color_vec, colors, idx, pntr, eps, nvars, nnz)
                            return(res)
                        } else {
                            stop("sparseHessianFD object not initialized")
                            res <- NULL
                        }
                        return(res)
                    },

                    fn = function(x) {
                        "Return function value, evaluated at x: fn(x, ...)"
                        fn1(x)
                    },

                    gr = function(x) {
                        "Return gradient, evaluated at x:  gr(x,...)"
                        gr1(x)
                    },

                    fngr = function(x) {
                        "Return list of function value and gradient, evaluated at x"
                        list(fn = fn(x), gr = gr(x))
                    },

                    fngrhs = function(x) {
                        "Return list of function value, gradient, and Hessian, evaluated at x"
                        list(fn = fn(x),
                             gr = gr(x),
                             hessian=hessian(x))
                    },


                    pointers = function(index1=TRUE) {
                        "Return list with indices (idx) and pointers (pntr) for sparsity pattern of the compressed sparse Hessian.  Since the Hessian is symmetric, the indices and pointers for row-oriented and column-oriented storage patterns are the same."
                        if (ready){
                            res <- list(idx = idx + as.integer(index1),
                                        pntr = pntr + as.integer(index1)
                                        )
                        } else {
                            stop("sparseHessianFD object not initialized")
                            res <- NULL
                        }
                        return(res)
                    },

                    coord2vec = function(j) {
                        z <- rep(0,nvars)
                        z[ind] <- eps
                        return(z)
                    })
                )


sparseHessianFD$usingMethods(coord2vec)
