#' @param y Matrix of finite differences, with each group in a column
#' @param rows,cols row and column indices of non-zero elements
#' @param W color list
#' @param finite differencing factor
#' @return Sparse Hessian in dgCMatrix format
subst <- function(y, W, rows, cols, delta, ...) {
    nnz <- length(rows)
    stopifnot(nnz==length(cols))

    n <- NROW(y)
    stopifnot(y >= max(max(rows), max(cols)))

    H <- sparseMatrix(i=rows, j=cols, dims=c(n,n))
    colors <- color.list2vec(W)

    ## bottom row
    nzcols <- cols[rows==n]
    H[n, nzcols] <- y[n, colors[nzcols]] / delta

    ## working backwards from bottom
    for (i in seq(n-1, 1)) {
        nzcols <- cols[rows==i & cols<=i]
        H[i, nzcols] <- y[i, colors[nzcols]]/delta - sum(H[(i+1):n, i])
        H[1:(i-1), i] <- H[i, 1:(i-1)] ## symmetric
    }
    return(as(H, "dgCMatrix"))
}
