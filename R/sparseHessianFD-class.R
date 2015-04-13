## This file is part of the sparseHessianFD package
## Copyright (C) 2015 Michael Braun
##
##' @title sparseHessianFD
##' @name sparseHessianFD-class
##' @description A reference class for computing sparse Hessians
##' @docType class
##' @field fn1 A closure for calling fn(x, ...).
##' @field gr1 A closure for calling gr(x, ...).
##' @field iRow,jCol Numeric vectors with row and column indices of the non-zero elements in the lower triangle (including diagonal) of the Hessian.
##' @field eps The perturbation amount for finite differencing of the gradient to compute the Hessian. Defaults to sqrt(.Machine$double.eps).
##' @field index1 TRUE if rows and cols use 1-based (R format) indexing (FALSE for 0-based (C format) indexing.
##' @field D raw finite differences (internal use only)
##' @field nvars Number of variables (length of x)
##' @field nnz Number of non-zero elements in the lower triangle of the Hessian.
##' @field ready TRUE if object has been initialized, and Hessian has been partitioned.
##' @field idx,pntr Column indices and row pointers for non-zero elements in lower triangle of the Hessian.  Row-oriented compressed storage.
##' @field colors A list representation of the partitioning of the columns of the Hessian.  This is used as part of the estimation algorithm.  Specifically, each list element is a "color," and each element in the corresponding vector is a column with that color.  See references.
##' @field colors_vec A vector representation of the partitioning of the columns.  There are nvars elements, one for each column of the Hessian.  The value corresponds to the "color" for that column.
##' @export sparseHessianFD
sparseHessianFD <-
    setRefClass("sparseHessianFD",
                fields=list(
                    fn1 = "function",
                    gr1 = "function",
                    iRow = "numeric",
                    jCol = "numeric",
                    eps = "numeric",
                    index1 = "logical",
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
                                          index1 = TRUE, ...) {
                        "Initialize object with functions to compute the objective function and gradient (fn and gr), row and column indices of non-zero elements (rows and cols), an initial variable vector x at which fn and gr can be evaluated, a finite differencing parameter eps, flags for 0 or 1-based indexing (index1), whether sparsity pattern is just for the lower triangle (indexLT), and other arguments (...) to be passed to fn and gr."

                        if (!is.null(direct)) {
                            warning(" 'direct' argument is ignored. Only indirect method is, and will be, supported.")
                        }
                        validate(fn, gr, rows, cols, x, eps, index1, ...)
                        ww <- which(cols <= rows)

                        initFields(fn1 = function(y) fn(y, ...),
                                   gr1 = function(y) gr(y, ...),
                                   iRow = as.integer(rows[ww]),
                                   jCol = as.integer(cols[ww]),
                                   eps = eps,
                                   index1 = index1,
                                   nvars = length(x),
                                   nnz = length(iRow),
                                   ready = FALSE)
                        if (any(cols > rows)) {
                            warning("Some elements are in upper triangle, and will be ignored.")
                        }



                        ptr <- Coord.to.Pointers(iRow, jCol, c(nvars, nvars),
                                                 triangle=TRUE,
                                                 lower=TRUE,
                                                 order="row",
                                                 index1=index1)

                        idx <<- ptr[[1]]
                        pntr <<- ptr[[2]]

                        ptr2 <- Coord.to.Pointers(iRow, jCol, c(nvars, nvars),
                                                 triangle=TRUE,
                                                 lower=TRUE,
                                                 order="symmetric",
                                                 index1=index1)

                        idx2 <- ptr2[[1]]
                        pntr2 <- ptr2[[2]]

                        colors <<- color_graph(pntr2-index1, idx2-index1, nvars)
                        colors_vec <<- rep(0L, nvars)
                        for (i in 1:length(colors)) {
                            colors_vec[colors[[i]]+1] <<- as.integer(i-1)
                        }

                        coord2vec <- function(j) {
                            z <- rep(0,nvars)
                            z[j+1] <- eps
                            return(z)
                        }


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
                            res <- NULL
                        }
                        return(res)
                    },

                    validate = function(fn, gr, rows, cols, x,
                                        eps, index1, ...) {

                        stopifnot(is.numeric(x),
                                  is.function(fn),
                                  is.function(gr),
                                  !is.null(rows),
                                  !is.null(cols),
                                  length(rows)==length(cols),
                                  all(is.finite(rows)),
                                  all(is.finite(cols)))

                        gradient <- gr(x,...)
                        val <- fn(x,...)

                        I1 <- as.integer(index1)
                        check.index1.row <- min(rows) >= I1 &
                          (max(rows) <= (nvars + I1 - 1))
                        check.index1.col <- min(cols) >= I1 &
                          (max(cols) <= (nvars + I1 - 1))

                        stopifnot(all(is.finite(gradient)),
                                  length(gradient)==nvars,
                                  check.index1.row,
                                  check.index1.col,
                                  is.finite(val)
                                  )
                    },

                    fd = function(d, x, grad.x) {
                        gr1(x+d) - grad.x
                    },

                    hessian = function(x) {
                        "Return sparse Hessian, evaluated at x, as a dgCMatrix object."
                        usingMethods(fd)
                        if (ready) {
                            grad.x <- gr1(x)
                            Y <- apply(D, 2, fd, x = x, grad.x = grad.x)
                            res <- subst(Y, colors_vec, colors, idx-index1, pntr-index1, eps, nvars, nnz)
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

                    pointers = function(out.index1=index1) {
                        "Return list with indices (idx) and pointers (pntr) for sparsity pattern of the compressed sparse Hessian.  Since the Hessian is symmetric, the indices and pointers for row-oriented and column-oriented storage patterns are the same."
                        if (ready){
                            res <- list(idx = idx + out.index1 - index1,
                                        pntr = pntr + out.index1 - index1
                                        )
                        } else {
                            stop("sparseHessianFD object not initialized")
                            res <- NULL
                        }
                        return(res)
                    })
                )


