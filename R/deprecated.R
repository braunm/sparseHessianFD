## deprecated.R -- Part of the sparseHessianFD package 
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.



## Place to hold functions that will not be maintained in the future


#' @name sparseHessianFD-deprecated
#' @aliases Sym.CSC.to.Matrix Coord.to.Sym.Pattern.Matrix
#' @title Deprecated functions
#' @description These functions were in earlier versions, but will no
#' longer be maintained, and are not even guaranteed to work now.
NULL

#' @title Sym.CSC.to.Matrix
#' @description Build sparse matrix from data in CSC (column
#' compressed) format.
#' @param H a list containing Hessian data.  See details.
#' @param nvars the number of rows (and columns) in the matrix.
#' @return An object of Matrix class.
#' @details Use Matrix::sparseMatrix instead of Sym.CSC.to.Matrix.
#' @rdname sparseHessianFD-deprecated
#' @export
Sym.CSC.to.Matrix <- function(H,nvars) {

    .Deprecated("Matrix::spMatrix")
    
  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)
}



#' @name Coord.to.Sym.Pattern.Matrix
#' @inheritParams Sym.CSC.to.Matrix
#' @details
#' Use Coord.to.Pattern.Matrix with symmetric=TRUE instead of Coord.to.Sym.Pattern.Matrix.
#' @rdname sparseHessianFD-deprecated
#' @export
Coord.to.Sym.Pattern.Matrix <- function(H, nvars) {

    .Deprecated("Coord.to.Pattern.Matrix")

## coords are for lower triangle, but coerces to symmetric pattern matrix
## H is a list with two integer vectors:  iRow and jCol

  
    
    res <- new("nsTMatrix",i=as.integer(H$iRow-1), j=as.integer(H$jCol-1),
               Dim=c(as.integer(nvars), as.integer(nvars)),uplo="L")

    return(res)

}

#' @name new.sparse.hessian.obj
#' @title Deprecated constructor
#' @param x variable vector for initialization
#' @param fn R function that returns function value
#' @param gr R function that returns the gradient of the function
#' @param hs list of two vectors:  row and column indices of non-zero
#' elements of lower triangle of Hessian.  See details.
#' @param fd.method If TRUE, use direct method for computatation.  Otherwise, use
#' indirect/substitution method.  See references.
#' @param eps The perturbation amount for finite differencing of the
#' gradient to compute the Hessian. Defaults to
#' sqrt(.Machine$double.eps).
#' @param ... Other parameters to be passed to fn and gr.
#' @details hs is a list of two elements:
#' \describe{
#' \item{iRow}{ Integer vector of row indices of non-zero elements in
#' lower triangle of Hessian.}
#' \item{jCol}{ Integer vector of column indices of non-zero elements in
#' lower triangle of Hessian.}
#' }
#' @rdname sparseHessianFD-deprecated
#' @export
new.sparse.hessian.obj <- function(x, fn, gr, hs, fd.method=0L, eps=sqrt(.Machine$double.eps),...) {

    .Deprecated("sparseHessianFD.new")
    if (is.null(hs))
        stop("sparseHessianFD: you must supply structure of the Hessian.")
    if (!is.list(hs))
        stop ("sparseHessianFD:  hs must be a list")
    if (!all.equal(names(hs), c("rows","cols"))) {
        if (all.equal(names(hs), c("iRow","jCol"))) {
            names(hs) <- c("rows","cols")
        }
    }
    if (!all.equal(names(hs),c("rows","cols"))) {
        stop ("sparseHessianFD:  Names of hs must be either (\"iRow, jCol\") or (\"rows, cols\")")
    }
    direct <- as.logical(fd.method)        
    return(sparseHessianFD.new(x, fn, gr, hs$rows, hs$cols, direct, eps, ...))
}
