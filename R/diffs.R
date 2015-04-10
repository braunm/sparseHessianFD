
#' @return vector of length n, with delta in each element
#' indexed by j, and zeros elsewhere.
#' @param j Indices of non-zero elements
#' @param n Length of vector
#' @param delta Value to be placed in each element indexed by j
.coord2vec <- function(j, n, delta=1) {
    ind <- as.integer(j)
    z <- rep(0,n)
    z[ind] <- delta
    return(z)
}

#' Compute finite difference
get.grad.delta <- function(d, x, df, gr, ...) {
    df(x + d, ...) - gr
}


## #' @param x Value at which to evaluate Hessian
## #' @param df function that returns gradient at x
## #' @param W coloring scheme from color.cols function
## #' @param delta finite differencing factor
## #' @return Finite difference gradients, with each group in a column.
## get.diffs <- function(x, df, W, delta, ...) {
##     k <- length(W) ## number of colors
##     n <- length(x) ## number of variables

##     gr <- df(x, ...) ## gradient at x

##     D <- sapply(W, .coord2vec, n, delta)

##     ## return gr(x+d) - gr(x) for each column of D
##     apply(D, 2, get.grad.delta, x=x, df=df, gr=gr)
## }


#' @param x Value at which to evaluate Hessian
#' @param df function that returns gradient at x
#' @param rows,cols row and column indices of non-zero elements in the lower
#' triangle of the sparse Hessian matrix.
#' @param W coloring scheme from color.cols function
#' @param delta finite differencing factor
#' @return Finite difference gradients, with each group in a column.
get.fd <- function(x, df, rows, cols, W, delta, ...) {

    nvars <- length(x)
    stopifnot(nvars >= max(max(rows), max(cols)))
    nnz <- length(rows)


    ## things that can happen at initialization
    M <- sparseMatrix(i=rows, j=cols)
    ptr <- Matrix.to.Pointers(M, order="row")
    colors <- as.integer(color.list2vec(W, nvars))

    gr <- df(x, ...) ## gradient at x
    D <- sapply(W, .coord2vec, nvars, delta)

    ## return gr(x+d) - gr(x) for each column of D
    Y <- apply(D, 2, get.grad.delta, x=x, df=df, gr=gr, ...)
    H <- subst2(Y, colors, W, ptr$jCol, ptr$ipntr, delta, nvars, nnz)


    return(H)
}
