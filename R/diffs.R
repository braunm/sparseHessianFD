
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


#' @description Returns row (column) indices for non-zero elements in
#' column (row) k.
#' @param k row or column for which you want the column or row indices of the
#' non-zero elements
#' @param a list representing the sparsity structure in column (row)-compressed
#' format. The first element contains the row (column) indices, and the second
#' contains the column (row) pointers.
#' @return An integer vector of row (column) indices.
#' @details All indices are 1-based (as in R).  The Matrix.to.Pointers function
#' will return a list suitable for use in this function.
get.indices <- function(k, Q) {
    Q[[1]][(Q[[2]][k]):(Q[[2]][k+1]-1)]
}


#' @param x Value at which to evaluate Hessian
#' @param df function that returns gradient at x
#' @param rows,cols row and column indices of non-zero elements in the lower
#' triangle of the sparse Hessian matrix.
#' @param W coloring scheme from color.cols function
#' @param delta finite differencing factor
#' @return Finite difference gradients, with each group in a column.
get.fd <- function(x, df, rows, cols, W, delta, ...) {

    k <- length(W)
    nvars <- length(x)
    stopifnot(nvars >= max(max(rows), max(cols)))
    nnz <- length(rows)


    ## things that can happen at initialization
    M <- sparseMatrix(i=rows, j=cols)
    ptr <- Matrix.to.Pointers(M, order="row")
    Sp <- lapply(1:nvars, get.indices, Q=ptr)
    colors <- as.integer(color.list2vec(W, nvars))
    colsize <- diff(ptr$ipntr)

    gr <- df(x, ...) ## gradient at x
    D <- sapply(W, .coord2vec, nvars, delta)

    ## return gr(x+d) - gr(x) for each column of D
    Y <- apply(D, 2, get.grad.delta, x=x, df=df, gr=gr)


    ##    tt <- system.time(H <- subst_C(Y, colors, W, Sp, colsize, delta))
    tt <- system.time(H <- subst2(Y, colors, W, Sp, delta, nnz))


    return(H)
}
