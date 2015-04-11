

context("colors")
test_that("colors", {
    library(dplyr)

    ## Test matrices H and L
    b <- 3
    n <- 200
    p <- 2

    H <- Matrix(FALSE, nrow=n*b+p, ncol=n*b+p)
    for (j in 1:b) {
        for (i in 1:n) {
            H[seq(i,n*b,by=n),(j-1)*n+i] <- TRUE
            H[(n*b+1):(n*b+p), (j-1)*n+i] <- TRUE
        }
    }
    H[,(n*b+1):(n*b+p)] <- T

    L <- bdiag(replicate(n, matrix(T,b,b), simplify=FALSE)) %>%
      rBind(matrix(TRUE, p,n*b)) %>%
      cBind(matrix(TRUE, b*n+p,p))

    testmat <- H
    hs <- Matrix.to.Coord(tril(testmat))


    nvars <- NROW(testmat)

##    W <- color.cols(hs$rows, hs$cols)
##    print("test_colors_R")
##    print(tril(testmat))
##    print(W)

    ##    print("test_colors_C")
    Z <- color.cols.C(hs$rows, hs$cols, nvars)
##    print(Z)


})
