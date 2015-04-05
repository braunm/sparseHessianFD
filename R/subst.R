#' @param y Matrix of finite differences, with each group in a column
#' @param rows,cols row and column indices of non-zero elements
#' @param W color list
#' @param finite differencing factor
#' @return Sparse Hessian in dgCMatrix format
subst <- function(y, W, rows, cols, delta, ...) {

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

        nzcols <- cols[rows==i & cols<=i]
        grp <- colors[i]
        yi <- y[nzcols, grp]
        ##     cat("i = ",i,"\n")
        ##     cat("\tnzcols = ",nzcols,"\n")
        ##     cat("\tgroups = ", grp, "\n")
        ##     cat("\tyi = ", yi, "\n")
        ind <- (seq(1,nvars) > i) & (colors==grp)
        H[i, nzcols] <- yi/delta - sum(H[ind, i])
        H[1:(i-1), i] <- H[i, 1:(i-1)] ## symmetric

    }
    return(H)
}
