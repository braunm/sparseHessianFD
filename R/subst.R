#' @param y Matrix of finite differences, with each group in a column
#' @param rows,cols row and column indices of non-zero elements
#' @param W color list
#' @param finite differencing factor
#' @return Sparse Hessian in dgCMatrix format
subst <- function(y, W, rows, cols, delta) {

    nnz <- length(rows)
    stopifnot(nnz==length(cols))

    nvars <- NROW(y)
    stopifnot(nvars >= max(max(rows), max(cols)))

    H <- as(sparseMatrix(i=rows, j=cols, dims=c(nvars, nvars)),
            "dgCMatrix")
    colors <- color.list2vec(W, nvars)
    ## bottom row
    nzcols <- cols[rows==nvars]
    yi <- y[, colors[nvars]]

    H[nvars, nzcols] <- yi/delta
    H[1:(nvars-1),nvars] <- H[nvars,1:(nvars-1)]

    ## working backwards from bottom
    for (i in seq(nvars-1, 1)) {
        grp <- colors[i]
        nzcols <- cols[rows==i & cols<=i]
        yi <- y[nzcols, grp]
        ind <- (seq(1,nvars) > i) & (colors==grp)
        H[i, nzcols] <- yi/delta - sum(H[ind, i])
        H[1:(i-1), i] <- H[i, 1:(i-1)] ## symmetric
    }
    return(H)
}

#' @param y Matrix of finite differences, with each group in a column
#' @param rows,cols row and column indices of non-zero elements
#' @param W color list
#' @param finite differencing factor
#' @return Sparse Hessian in dgCMatrix format
subst.C <- function(Y, W, rows, cols, delta) {

    nnz <- length(rows)
    stopifnot(nnz==length(cols))

    nvars <- NROW(Y)
    stopifnot(nvars >= max(max(rows), max(cols)))

    colors <- as.integer(color.list2vec(W, nvars))
    Sp <- vector("list", length=nvars)

    for (i in 1:nvars) {
        q <- which(rows==i)
        Sp[[i]] <- as.integer(sort(cols[q]))
    }

    H <- subst2(Y, colors, W, Sp, delta)
    return(H)

}
