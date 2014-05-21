## matrices.R --   Part of the sparseHessianFD package for the R programming language.
##
## Copyright (C) 2013 Michael Braun



Sym.CSC.to.Matrix <- function(H,nvars) {

## H is a list of data stored in compressed sparse column (CSC) format
## returns a sparse Matrix object

  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)

}


Matrix.to.Coord <- function(M) {

  res <- vector("list",length=2)
  names(res) <- c("iRow","jCol")
  M <- tril(as(M,"TsparseMatrix"))
  res$iRow <- as.integer(M@i)+1 ## return to 1-based indexingk
  res$jCol <- as.integer(M@j)+1
  return(res)

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
