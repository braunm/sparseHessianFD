## Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

#' @name coloring
#' @title Triangular partitioning of variables
#' @param L  nMatrix with sparsity pattern of lower triangle of Hessian
#' @return Integer vector of length nvars with color assignments for each variable.
#' @description cyclic coloring from a lower triangular pattern matrix
#' @export
coloring <- function(L) {

    stopifnot(is(L,"nMatrix"),
              Matrix::isTriangular(L, upper=FALSE)
              )

    nvars <- NROW(L)
    L <- as(L,"nMatrix")
    G <- Matrix::crossprod(L)  # intersection graph
 ##   ptr <- Matrix.to.Pointers(G, order="symmetric", index1=FALSE)

    ptr <- Matrix.to.Pointers(G, index1=FALSE)

    colors_vec <- get_colors(ptr[[2]], ptr[[1]], nvars)
    return(colors_vec)
}



