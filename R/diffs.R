
#' @return vector of length n, with delta in each element
#' indexed by j, and zeros elsewhere.
#' @param j Indices of non-zero elements
#' @param n Length of vector
#' @param delta Value to be placed in each element indexed by j
.coord2vec <- function(j, n, delta=1) {
    ind <- as.integer(j)
    z <- rep(0,n)
    z[ind] <- 1
    return(z)
}

#' Compute finite difference
get.grad.delta <- function(d, x, df, gr, ...) {
    df(x + d, ...) - gr
}


#' @param x Value at which to evaluate Hessian
#' @param df function that returns gradient at x
#' @param rows row indices of non-zero elements
#' @param cols col indices of on-zero elements
#' @param W coloring scheme from color.cols function
#' @param delta finite differencing factor
#' @return Finite difference gradients, with each group in a column.
get.diffs <- function(x, df, rows, cols, W, delta, ...) {
    k <- length(W) ## number of colors
    n <- length(x) ## number of variables
    D <- sapply(W, .coord2vec, n, delta)

    gr <- df(x, ...) ## gradient at x

    ## return gr(x+d) - gr(x) for each column of D
    apply(D, 2, get.grad.delta, x=x, df=df, gr=gr)
}
