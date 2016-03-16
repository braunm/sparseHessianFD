

make_pointers <- function(rows, cols, nvars, indexLT=TRUE, index1=TRUE, out.index1=index1) {

    stopifnot(is.logical(indexLT),
              is.logical(index1),
              is.logical(out.index1)
              )

    if (indexLT) {
        R <- as(sparseMatrix(i=rows, j=cols, dims=c(nvars, nvars), index1=index1),
                "nMatrix") ## LT
        C <- as(sparseMatrix(i=cols, j=rows, dims=c(nvars, nvars), index1=index1),
                "nMatrix") ## UT
        A <- as(R + C, "ngCMatrix") ## symmetric , but stored as general CSC sparse
    } else {
        A <- as(sparseMatrix(i=rows, j=cols, dims=c(nvars, nvars)), "ngCMatrix")
        if (!isSymmetric(A)) {
            stop("indexLT == FALSE, but matrix is not symmetric")
        }
    }

    list(idx = A@i + as.integer(out.index1),
         pntr = A@p + as.integer(out.index1))

}


