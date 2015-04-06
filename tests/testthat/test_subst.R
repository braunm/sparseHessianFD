

context("subst_C")
test_that("subst_C", {

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

    Y1 <- get.diffs(P, df=f1$gr, pattern1$rows, pattern1$cols, W1, delta)
    Y2 <- get.diffs(P, df=f2$gr, pattern2$rows, pattern2$cols, W2, delta)


    H1 <- subst.C(Y1, W1, pattern1$rows, pattern1$cols, delta)
    H2 <- subst.C(Y2, W2, pattern2$rows, pattern2$cols, delta)

browser()

    expect_equal(true.hess1, H1, tolerance=1e-7)
    expect_equal(true.hess2, H2, tolerance=1e-7)

})
