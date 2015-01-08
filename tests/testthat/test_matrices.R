context("Matrix utility functions")

test_that("Matrix.to.Coord", {

    nnz <- 6
    k <- 5
    rows <- c(1:5,5)
    cols <- c(1:5,2)
    vals <- c(1,2,3,4,5,6)


    
    ## M <- Matrix(0,k,k)
    ## for (i in 1:nnz) {
    ##     M[rows[i], cols[i]] <- vals[i]
    ## }
    

    res <- Matrix.to.Coord(M)

    print(cbind(res$iRow, res$jCol))
    expect_equal(names(res), c("iRow","jCol"))
 ##   expect_equal(res$iRow, rows)
 ##   expect_equal(res$jCol, cols)
})
