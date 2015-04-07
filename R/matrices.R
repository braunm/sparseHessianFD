## matrices.R -- Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

#' @name Matrix.to.Coord
#' @title Row and column indices from sparse matrix.
#' @description Returns list of row and column indices of the non-zero
#' elements of a sparse matrix.
#' @param M A sparse Matrix, as defined in the Matrix package.
#' @return A list with two named elements.
#' \describe{
#' \item{rows}{ Integer vector containing row indices of non-zero elements}
#' \item{cols}{ Integer vector containing column indices of non-zero elements}
#' }
#' @export
Matrix.to.Coord <- function(M) {
  res <- vector("list",length=2)
  names(res) <- c("rows","cols")
  M <- as(M,"TsparseMatrix")
  res$rows <- as.integer(M@i + 1) ## return to 1-based indexing
  res$cols <- as.integer(M@j + 1)
  return(res)
}

#' @name Coord.to.Pattern.Matrix
#' @aliases Coord.to.Pattern.Matrix
#' @title Pattern matrix from row and column indices.
#' @description Converts row and column indices to a pattern Matrix
#' object of Matrix class
#' @param rows,cols row and column indices of non-zero elements
#' @param dims 2-element vector for number of rows and columns in
#' matrix
#' @param compressed If TRUE, returns a matrix is compressed column (default=TRUE)
#' @param symmetric If TRUE, matrix will be symmetric, and only the
#' lower triangular elements need to be provided (default=FALSE)
#' @return A sparse pattern matrix
#' @details This function is useful to prototype a sparsity pattern.
#' No assumptions are made about symmetry.
#' @export
Coord.to.Pattern.Matrix <- function(rows, cols, dims, compressed=TRUE,
                                    symmetric=FALSE) {

    res <- sparseMatrix(i=as.integer(rows),
                        j=as.integer(cols),
                        dims=dims,
                        giveCsparse=compressed,
                        symmetric=symmetric)
    return(res)

}


#' @name Matrix.to.Pointers
#' @title Column (row) indices and row (column) pointers from sparse matrix.
#' @description Returns list of column (row) indices row (column) pointers
#' of the non-zero elements of a sparse matrix stored in row (column)-oriented
#' compressed format.
#' @param M A sparse Matrix, as defined in the Matrix package.
#' @return A list with two named elements.
#' \describe{
#' \item{rows}{ Integer vector containing row indices of non-zero elements}
#' \item{cols}{ Integer vector containing column indices of non-zero elements}
#' }
#' @export
Matrix.to.Pointers <- function(M, order=c("column", "row")) {
    res <- vector("list",length=2)
    if (order == "row") {
        M <- as(M,"RsparseMatrix")
        names(res) <- c("jCol","ipntr")
        res$jCol <- as.integer(M@j + 1)
        res$ipntr <- as.integer(M@p + 1)
    } else {
        M <- as(M,"CsparseMatrix")
        names(res) <- c("iRow","jpntr")
        res$iRow <- as.integer(M@i + 1)
        res$jpntr <- as.integer(M@p + 1)
    }
    return(res)
}
