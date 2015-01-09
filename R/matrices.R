## matrices.R --   Part of the sparseHessianFD package for the R programming language.
##
## Copyright (C) 2013 Michael Braun

#' @name Matrix.to.Coord
#' @title Row and column indices from sparse symmetric matrix.
#' @description Returns list of row and column indices of the non-zero
#' elements of the lower triangle of a symmetric matrix.
#' @param M A symmetric Matrix
#' @details Takes an object of class Matrix and returns the row and
#' column indices of the non-zero elements in the lower triangle 
#' @return A list with two named elements.
#' \describe{
#' \item{iRow}{ Integer vector containing row indices of non-zero elements in the lower triangle of M}
#' \item{jCol}{ Integer vector containing column indices of non-zero
#' elements in the lower triangle of M}
#' }
#' @export
Matrix.to.Coord <- function(M) {
  res <- vector("list",length=2)
  names(res) <- c("iRow","jCol")
  M <- tril(as(M,"TsparseMatrix"))
  res$iRow <- as.integer(M@i + 1) ## return to 1-based indexing
  res$jCol <- as.integer(M@j + 1)
  return(res)
}


#' @title Sym.CSC.to.Matrix
#' @description Build sparse matrix from data in CSC (column
#' compressed) format.
#' @param H a list containing Hessian data.  See details.
#' @param nvars the number of rows (and columns) in the matrix.
#' @return An object of Matrix class.
#' @details
#' H is a list with three elements.
#' \describe{
#' \item{iRow}{Row indices for each of the non-zero elements in the lower triangle of H}
#' \item{jpntr}{A vector of length \code{nvars+1}.  jpntr[j] is the index of the element in vals that is the first non-zero element in the jth column of the matrix.}
#' \item{vals}{The values of the non-zero elements in the lower triangle of the matrix.}
#' }
#' The input list contains information about only the non-zero
#' elements in the lower triangle of the matrix.
#' @export
Sym.CSC.to.Matrix <- function(H,nvars) {
  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)
}

#' @name Coord.to.Pattern.Matrix
#' @aliases Coord.to.Pattern.Matrix
#' @title Pattern matrix from row and column indices.
#' @description Converts row and column indices to a pattern Matrix object of Matrix class
#' @param H a list containing matrix structure data.  See details.
#' @param nrows the number of rows in the matrix.
#' @param ncols the number of columns in the matrix.
#' @return An object of Matrix class.
#' @details H is a list with two vectors.
#' \describe{
#' \item{iRow}{Row indices for each of the non-zero elements in the matrix.}
#' \item{jCol}{Column indices for each of the non-zero elements in the matrix.}}
#' @export
Coord.to.Pattern.Matrix <- function(H,nrows, ncols=nrows) {

  ## H is a list with two integer vectors:  iRow and jCol
  
  res <- sparseMatrix(i=H$iRow,j=H$jCol, dims=c(as.integer(nrows), as.integer(ncols)))
  return(res)

}

Coord.to.Sym.Pattern.Matrix <- function(H,nvars) {

## coords are for lower triangle, but coerces to symmetric pattern matrix
## H is a list with two integer vectors:  iRow and jCol

  
  res <- new("nsTMatrix",i=as.integer(H$iRow-1), j=as.integer(H$jCol-1),
             Dim=c(as.integer(nvars), as.integer(nvars)),uplo="L")
  return(res)

}
