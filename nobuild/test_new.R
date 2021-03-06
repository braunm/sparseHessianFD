

context("newmethod")
test_that("newmethod", {

    set.seed(123)
    data(binary_large)
    binary <- binary_large

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
    pattern1LT <- Matrix.to.Coord(tril(true.hess1))
    pattern2LT <- Matrix.to.Coord(tril(true.hess2))

    pattern1 <- Matrix.to.Coord(true.hess1)
    pattern2 <- Matrix.to.Coord(true.hess2)

    W1 <- color.cols.C(pattern1LT$rows, pattern1LT$cols, nvars)
    W2 <- color.cols.C(pattern2LT$rows, pattern2LT$cols, nvars)

    delta <- 1e-7
##browser()

    H1 <- get.fd(P, df=f1$gr, pattern1LT$rows, pattern1LT$cols, W1, delta)
    H2 <- get.fd(P, df=f2$gr, pattern2LT$rows, pattern2LT$cols, W2, delta)

    expect_equal(true.hess1, H1)
    expect_equal(true.hess2, H2)

})
