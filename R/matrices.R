## matrices.R -- Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

#' @name Matrix.to.Coord
#' @title Row and column indices from sparse matrix.
#' @description Returns list of row and column indices of the non-zero
#' elements of a sparse matrix.
#' @param M A matrix that can be coerced into a \code{TsparseMatrix} object, as defined in the Matrix package.  Base R matrices, and most Matrix matrices, are valid.
#' @param index1 TRUE if the index of the first element should be 1, and FALSE if 0.
#' @return A list with two named elements.
#' \describe{
#' \item{rows}{ Integer vector containing row indices of non-zero elements}
#' \item{cols}{ Integer vector containing column indices of non-zero elements}
#' }
#' @export
Matrix.to.Coord <- function(M, index1=TRUE) {
  res <- vector("list",length=2)
  names(res) <- c("rows","cols")
  M <- as(M,"TsparseMatrix")
  res$rows <- as.integer(M@i) + as.integer(index1)
  res$cols <- as.integer(M@j) + as.integer(index1)
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
#' @param index1 TRUE if input row and col use 1-based indexing, and FALSE for 0-based indexing.
#' @return A sparse pattern matrix
#' @details This function is useful to prototype a sparsity pattern.
#' No assumptions are made about symmetry.
#' @export
Coord.to.Pattern.Matrix <- function(rows, cols, dims, compressed=TRUE,
                                    symmetric=FALSE, index1=TRUE) {

    res <- sparseMatrix(i=as.integer(rows),
                        j=as.integer(cols),
                        dims=dims,
                        giveCsparse=compressed,
                        symmetric=symmetric,
                        index1=index1)
    return(res)
}


#' @name Matrix.to.Pointers
#' @title Column (row) indices and row (column) pointers from sparse matrix.
#' @description Returns list of column (row) indices row (column) pointers
#' of the non-zero elements of a sparse matrix stored in row (column)-oriented
#' compressed format.
#' @param M A sparse Matrix, as defined in the Matrix package.
#' @param order Pointers can be ordered according to a row-oriented or column-oriented compression pattern.  If symmetric is chosen, the compression method does not matter.  The input matrix must be symmetric (although not necessarily stored that way) to get symmetric output.  The only difference between "symmetric," and either "row" or "column," is the names of the elements in the list.
#' @param index1 TRUE if return values use 1-based indexing, and FALSE for 0-based indexing.
#' @return A list with two named elements. If matrix is row-oriented,
#' \describe{
#' \item{jCol}{ Integer vector containing column indices of non-zero elements}
#' \item{ipntr}{ Integer vector containing pointers to elements of jCol at which the next row begins.}
#' }
#'  If matrix is column-oriented,
#'  \describe{
#' \item{iRow}{ Integer vector containing row indices of non-zero elements}
#' \item{jpntr}{ Integer vector containing pointers to elements of iRow at which the next column begins.}
#' }
#' If the matrix is symmetric, then the row/column orientation does not matter.
#'  \describe{
#' \item{idx}{ Integer vector containing indices of non-zero elements}
#' \item{pntr}{ Integer vector containing pointers to elements of idx at which the next row or column begins.}
#' }
#' @export
Matrix.to.Pointers <- function(M, order="symmetric",
                               index1=TRUE) {


    stopifnot(order %in% c("column","row","symmetric"))
    res <- vector("list",length=2)
    if (order == "symmetric") {
        stopifnot(Matrix::isSymmetric(M))
        M <- as(M,"dgCMatrix")
        names(res) <- c("idx","pntr")
        res$idx <- as.integer(M@i) + as.integer(index1)
        res$pntr <- as.integer(M@p) + as.integer(index1)
    } else {
        if (order == "row") {
            M <- as(M,"RsparseMatrix")
            names(res) <- c("jCol","ipntr")
            res$jCol <- as.integer(M@j) + as.integer(index1)
            res$ipntr <- as.integer(M@p) + as.integer(index1)
        } else {
            M <- as(M,"CsparseMatrix")
            names(res) <- c("iRow","jpntr")
            res$iRow <- as.integer(M@i) + as.integer(index1)
            res$jpntr <- as.integer(M@p) + as.integer(index1)
        }
    }
    return(res)
}


#' @name Coord.to.Pointers
#' @title Column (row) indices and row (column) pointers from sparse matrix.
#' @description Returns list of column (row) indices row (column) pointers
#' of the non-zero elements of a sparse matrix stored in row (column)-oriented
#' compressed format.
#' @param rows,cols row and column indices of non-zero elements
#' @param dims 2-element vector for number of rows and columns in
#' matrix
#' @param order Pointers can be ordered according to a row-oriented or column-oriented compression pattern.  If symmetric is chosen, the compression method does not matter.  The input matrix must be symmetric (although not necessarily stored that way) to get symmetric output.  For symmetric input, the only difference between "symmetric," and either "row" or "column," is the names of the elements in the list. See details.
#' @param triangle Is input intended to be a triangular (TRUE) or full (FALSE) matrix.  See Details for how this argument is interpreted for different values of \code{order}.
#' @param lower If \code{triangular} is true, this argument identifies if the input matrix is lower- or upper-triangular.  This argument is ignored if \code{triangle} is FALSE.
#' @param index1 TRUE if using 1-based indexing.  FALSE for 0-based indexing.
#' @details \code{triangle} and \code{order} have the following interpretation:
#' \describe{
#' \item{triangle=TRUE}{Input \code{rows} and {cols} represent lower or upper triangle of a matrix. If \code{order="symmetric"}, then the output list will be for a full, symmetric matrix. Otherwise, the output list will be for only the lower or upper triangle.  Any elements outside of the specified triangle will trigger an error.}
#' \item{triangle=FALSE}{Input \code{rows} and {cols} represent a full matrix. If that matrix is not symmetric, then \code{order=="symmetric"} will trigger an error.}
#' }
#' @return A list with two named elements.  See 'Matrix.to.Pointers'.
#' @export
Coord.to.Pointers <- function(rows, cols, dims,
                              triangle=TRUE, lower=TRUE,
                              order=c("column", "row", "symmetric"),
                              index1=TRUE) {

    stopifnot(is.logical(triangle),
              is.logical(index1),
              is.logical(lower),
              length(rows) == length(cols)
              )

    if (triangle) {
        if (lower) {
            stopifnot(all(rows >= cols))
        } else {
            stopifnot(all(rows <= cols))
        }
    }

    R <- as(sparseMatrix(i=rows, j=cols, dims=dims,
                         index1=index1, giveCsparse=TRUE), "nMatrix")
    if (triangle & order=="symmetric") {
        C <- as(sparseMatrix(i=cols, j=rows, dims=dims, index1=index1), "nMatrix")
        A <- as(R + C, "ngCMatrix") ## symmetric , but stored as general CSC sparse
    } else {
        A <- R
    }
    Matrix.to.Pointers(A, order, index1, index1)
}


