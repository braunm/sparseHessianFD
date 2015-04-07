

context("compare")
test_that("compare", {

    library(dplyr)
    set.seed(123)
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

    true.val1 <- f1$fn(P)
    true.grad1 <- f1$gr(P)
    true.hess1 <- drop0(f1$hessian(P))

    true.val2 <- f2$fn(P)
    true.grad2 <- f2$gr(P)
    true.hess2 <- drop0(f2$hessian(P))

    ## Get hessian structure
    pattern1 <- Matrix.to.Coord(tril(true.hess1))
    pattern2 <- Matrix.to.Coord(tril(true.hess2))

    W1 <- color.cols(pattern1$rows, pattern1$cols)
    W2 <- color.cols(pattern2$rows, pattern2$cols)


    delta <- 1e-7

    H1 <- get.fd(P, df=f1$gr, pattern1$rows, pattern1$cols, W1, delta)
    H2 <- get.fd(P, df=f2$gr, pattern2$rows, pattern2$cols, W2, delta)

    expect_equivalent(true.hess1, drop0(H1))
    expect_equivalent(true.hess2, drop0(H2))

    obj10 <- new("sparseHessianFD", nvars, f1$fn, f1$gr)
    obj10$hessian.init(pattern1$rows, pattern1$cols, 0, delta)

    obj11 <- new("sparseHessianFD", nvars, f1$fn, f1$gr)
    obj11$hessian.init(pattern1$rows, pattern1$cols, 1, delta)

    obj20 <- new("sparseHessianFD", nvars, f2$fn, f2$gr)
    obj20$hessian.init(pattern2$rows, pattern2$cols, 0, delta)

    obj21 <- new("sparseHessianFD", nvars, f2$fn, f2$gr)
    obj21$hessian.init(pattern2$rows, pattern2$cols, 1, delta)

    test.val10 <- obj10$fn(P)
    test.grad10 <- obj10$gr(P)
    test.hess10 <- obj10$hessian(P)

    test.val11 <- obj11$fn(P)
    test.grad11 <- obj11$gr(P)
    test.hess11 <- obj11$hessian(P)

    test.val20 <- obj20$fn(P)
    test.grad20 <- obj20$gr(P)
    test.hess20 <- obj20$hessian(P)

    test.val21 <- obj21$fn(P)
    test.grad21 <- obj21$gr(P)
    test.hess21 <- obj21$hessian(P)

    expect_equivalent(drop0(H1), test.hess10)
    expect_equivalent(drop0(H1), test.hess11)
    expect_equivalent(drop0(H2), test.hess20)
    expect_equivalent(drop0(H2), test.hess21)




})
