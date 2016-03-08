## Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

context("colors")

test_that("colors", {

    set.seed(123)
    data(binary_small)
    binary <- binary_small

    N <- length(binary$Y)
    k <- NROW(binary$X)
    nvars <- as.integer(N*k + k)
    P <- rnorm(nvars) ## random evaluation values

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

    g1 <- coloring(as(tril(true.hess1),"nMatrix"))
    g2 <- coloring(as(tril(true.hess2), "nMatrix"))

    pat1 <- Matrix.to.Coord(tril(true.hess1))
    pat2 <- Matrix.to.Coord(tril(true.hess2))

    obj1 <- sparseHessianFD(P, f1$fn, f1$gr, pat1$rows, pat1$cols)
    obj2 <- sparseHessianFD(P, f2$fn, f2$gr, pat2$rows, pat2$cols)


})





