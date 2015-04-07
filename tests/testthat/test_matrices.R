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



test_that("Matrix.to.Pointers", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,2)
    cols <- c(1,2,2,3,4,5,5)
    vals <- c(1,2,6,3,4,5,6)


    M1 <- sparseMatrix(i=rows, j=cols, x=vals, dims=c(k,k) )
    C1 <- Matrix.to.Pointers(M1, order="column")
    expect_equal(names(C1), c("iRow","jpntr"))

    M2 <- sparseMatrix(i=C1$iRow, p=C1$jpntr-1, dims=c(k,k))
    M3 <- as(M1, "nMatrix")
    expect_equal(M2, M3)

    M8 <- as(M1, "RsparseMatrix")
    M9 <- as(M8, "nMatrix")
    C4 <- Matrix.to.Pointers(M8, order="row")
    expect_equal(names(C4), c("jCol","ipntr"))
    M4 <- sparseMatrix(j=C4$jCol, p=C4$ipntr-1, dims=c(k,k))
    expect_equal(M4, M9)
})
