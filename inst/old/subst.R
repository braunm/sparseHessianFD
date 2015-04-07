## #' @param y Matrix of finite differences, with each group in a column
## #' @param rows,cols row and column indices of non-zero elements
## #' @param W color list
## #' @param finite differencing factor
## #' @return Sparse Hessian in dgCMatrix format
## subst_R <- function(y, W, rows, cols, delta) {

##     nnz <- length(rows)
##     stopifnot(nnz==length(cols))

##     nvars <- NROW(y)
##     stopifnot(nvars >= max(max(rows), max(cols)))

##     H <- as(sparseMatrix(i=rows, j=cols, dims=c(nvars, nvars)),
##             "dgCMatrix")
##     colors <- color.list2vec(W, nvars)
##     ## bottom row
##     nzcols <- cols[rows==nvars]
##     yi <- y[, colors[nvars]]

##     H[nvars, nzcols] <- yi/delta
##     H[1:(nvars-1),nvars] <- H[nvars,1:(nvars-1)]

##     ## working backwards from bottom
##     for (i in seq(nvars-1, 1)) {
##         grp <- colors[i]
##         nzcols <- cols[rows==i & cols<=i]
##         yi <- y[nzcols, grp]
##         ind <- (seq(1,nvars) > i) & (colors==grp)
##         H[i, nzcols] <- yi/delta - sum(H[ind, i])
##         H[1:(i-1), i] <- H[i, 1:(i-1)] ## symmetric
##     }
##     return(H)
## }

## #' @param y Matrix of finite differences, with each group in a column
## #' @param rows,cols row and column indices of non-zero elements
## #' @param W color list
## #' @param finite differencing factor
## #' @return Sparse Hessian in dgCMatrix format
## subst <- function(Y, W, rows, cols, delta) {

##     nnz <- length(rows)
##     stopifnot(nnz==length(cols))

##     nvars <- NROW(Y)
##     stopifnot(nvars >= max(max(rows), max(cols)))

##     colors <- as.integer(color.list2vec(W, nvars))

##     ## Cols do not need to be sorted
##     Sp <- lapply(1:nvars, function(i) as.integer(cols[rows==i]))
##   ##  Sp <- lapply(1:nvars, function(i) sort(as.integer(cols[rows==i])))

##     H <- subst_C(Y, colors, W, Sp, delta)
##     return(H)

## }
