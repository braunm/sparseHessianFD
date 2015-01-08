## matrices.R --   Part of the sparseHessianFD package for the R programming language.
##
## Copyright (C) 2013 Michael Braun



#' @name Matrix.to.Coord
#' @title Row and column indices of non-zero elements of a sparse symmetric matrix
#' @param M A symmwetric Matrix
#' @description Takes an object of class Matrix and returns the row and column indices of the non-zero elements in the lower triangle
#' @return A list with two named elements.
#' \item{iRow} Integer vector containing row indices of non-zero elements in the lower triangle of M
#' \item{jCol} Integer vector containing column indices of non-zero elements in the lower triangle of M
#' @export
Matrix.to.Coord <- function(M) {
  res <- vector("list",length=2)
  names(res) <- c("iRow","jCol")
  M <- tril(as(M,"TsparseMatrix"))
  res$iRow <- as.integer(M@i + 1) ## return to 1-based indexingk
  res$jCol <- as.integer(M@j + 1)
  return(res)
}

Sym.CSC.to.Matrix <- function(H,nvars) {

## H is a list of data stored in compressed sparse column (CSC) format
## returns a sparse Matrix object

  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)

}




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
