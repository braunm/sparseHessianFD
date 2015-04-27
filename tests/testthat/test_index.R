## Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

context("indexing")

test_that("indexing", {

    set.seed(1234)
    data(binary_small)
    binary <- binary_small

    N <- length(binary$Y)
    k <- NROW(binary$X)
    nvars <- as.integer(N*k + k)
    P <- rnorm(nvars) ## random starting values

    Omega <- diag(k)
    priors <- list(inv.Omega = solve(Omega),
                   inv.Sigma = rWishart(1,k+5,diag(k))[,,1])

    make.funcs <- function(D, priors, order.row) {
        res <- vector("list", length=3)
        names(res) <- c("fn", "gr", "hessian")
        res$fn <-  function(pars) {
            binary.f(pars, data=D, priors=priors, order.row=order.row)
        }
        res$gr <-  function(pars) {
            binary.grad(pars, data=D, priors=priors, order.row=order.row)
        }
        res$hessian <-  function(pars) {
            binary.hess(pars, data=D, priors=priors, order.row=order.row)
        }
        return(res)
    }

    f1 <- make.funcs(D=binary, priors=priors, order.row=FALSE) ## block-arrow
    f2 <- make.funcs(D=binary, priors=priors, order.row=TRUE) ## off-diagonals

    ## True values for test


    true.hess1 <- drop0(f1$hessian(P))
    true.hess2 <- drop0(f2$hessian(P))

    ## Get hessian structure
    pat1L <- Matrix.to.Coord(tril(true.hess1))
    pat2L <- Matrix.to.Coord(tril(true.hess2))
    pat1S <- Matrix.to.Coord(true.hess1)
    pat2S <- Matrix.to.Coord(true.hess2)



    obj1L <- sparseHessianFD(P, f1$fn, f1$gr, pat1L$rows, pat1L$cols,
                index1=TRUE)
    obj2L <- sparseHessianFD(P, f2$fn, f2$gr, pat2L$rows, pat2L$cols,
                index1=TRUE)
    obj3L <- sparseHessianFD(P, f1$fn, f1$gr, pat1L$rows-1, pat1L$cols-1,
                index1=FALSE)
    obj4L <- sparseHessianFD(P, f2$fn, f2$gr, pat2L$rows-1, pat2L$cols-1,
                 index1=FALSE)
    ## obj1S <- sparseHessianFD(P, f1$fn, f1$gr, pat1S$rows, pat1S$cols,
    ##              index1=TRUE)
    ## obj2S <- sparseHessianFD(P, f2$fn, f2$gr, pat2S$rows, pat2S$cols,
    ##              index1=TRUE)
    ## obj3S <- sparseHessianFD(P, f1$fn, f1$gr, pat1S$rows-1, pat1S$cols-1,
    ##              index1=FALSE)
    ## obj4S <- sparseHessianFD(P, f2$fn, f2$gr, pat2S$rows-1, pat2S$cols-1,
    ##              index1=FALSE)

    print(true.hess1[c(9,10,1:8),c(9,10,1:8)])
    H1L <- obj1L$hessian(P)



    H2L <- obj2L$hessian(P)
    H3L <- obj3L$hessian(P)
    H4L <- obj4L$hessian(P)
    ## H1S <- obj1S$hessian(P)
    ## H2S <- obj2S$hessian(P)
    ## H3S <- obj3S$hessian(P)
    ## H4S <- obj4S$hessian(P)


    expect_equal(H1L, true.hess1, tolerance=5e-8)
    expect_equal(H2L, true.hess2, tolerance=5e-8)
    expect_equal(H1L, H3L)
    expect_equal(H2L, H4L)
    ## expect_equal(H1S, true.hess1, tolerance=5e-8)
    ## expect_equal(H2S, true.hess2, tolerance=5e-8)
    ## expect_equal(H1S, H3S)
    ## expect_equal(H2S, H4S)
})


