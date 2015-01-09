## test_matrices.R -- Part of the sparseHessianFD package 
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.


context("Matrix utility functions")

test_that("Matrix.to.Coord", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,2)
    cols <- c(1,2,2,3,4,5,5)
    vals <- c(1,2,6,3,4,5,6)


    M1 <- sparseMatrix(i=rows, j=cols, x=vals, dims=c(k,k) )
    C1 <- Matrix.to.Coord(M1)
    expect_equal(names(C1), c("rows","cols"))
    M2 <- sparseMatrix(i=C1$rows, j=C1$cols, dims=c(k,k))
    M3 <- as(M1, "nMatrix")
    expect_equal(M2, M3)

 
    P1 <- Coord.to.Pattern.Matrix(C1$rows, C1$cols, c(k,k))
    expect_is(P1, "ngCMatrix")
    T1 <- tril(M3)
    expect_equal(P1, M3)
    expect_equal(T1, tril(M3))

    C2 <- Matrix.to.Coord(tril(M1))
    P2 <- Coord.to.Pattern.Matrix(C2$rows, C2$cols, c(k,k), symmetric=TRUE)
    expect_that(P2, is_a("nsCMatrix"))
    expect_equal(forceSymmetric(tril(P1)), P2)

    P3 <- Coord.to.Pattern.Matrix(C2$rows, C2$cols, c(k,k), compressed=FALSE)
    expect_is(P3, "ngTMatrix")
    expect_equivalent(P1, P3)

})
