library(Matrix)
library(testthat)

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
true.hess1 <- drop0(f1$hessian(P))
true.hess2 <- drop0(f2$hessian(P))

## Get hessian structure
pattern1 <- Matrix.to.Coord(tril(true.hess1))
pattern2 <- Matrix.to.Coord(tril(true.hess2))

time.col1 <- system.time(W1 <- color.cols(pattern1$rows, pattern1$cols))
time.col2 <- system.time(W2 <- color.cols(pattern2$rows, pattern2$cols))

delta <- 1e-7

time.diff1 <- system.time(Y1 <- get.diffs(P, df=f1$gr, pattern1$rows, pattern1$cols, W1, delta))
time.diff2 <- system.time(Y2 <- get.diffs(P, df=f2$gr, pattern2$rows, pattern2$cols, W2, delta))

time.subst1 <- system.time(H1 <- subst(Y1, W1, pattern1$rows, pattern1$cols, delta))
time.subst2 <- system.time(H2 <- subst(Y2, W2, pattern2$rows, pattern2$cols, delta))

time.set.10 <- system.time(obj10 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))
time.init.10 <- system.time(obj10$hessian.init(pattern1$rows, pattern1$cols, 0, delta))

time.set.11 <- system.time(obj11 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))
time.init.11 <- system.time(obj11$hessian.init(pattern1$rows, pattern1$cols, 1, delta))

time.set.20 <- system.time(obj20 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))
time.init.20 <- system.time(obj20$hessian.init(pattern2$rows, pattern2$cols, 0, delta))

time.set.21 <- system.time(obj21 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))
time.init.21 <- system.time(obj21$hessian.init(pattern2$rows, pattern2$cols, 1, delta))


time.fd.10 <- system.time(test.hess10 <- obj10$hessian(P))
time.fd.11 <- system.time(test.hess11 <- obj11$hessian(P))
time.fd.20 <- system.time(test.hess20 <- obj20$hessian(P))
time.fd.21 <- system.time(test.hess21 <- obj21$hessian(P))


expect_equal(true.hess1, H1, tolerance=1e-6)
expect_equal(true.hess2, H2, tolerance=1e-6)
expect_equal(H1, test.hess10)
expect_equal(H1, test.hess11)
expect_equal(H2, test.hess20)
expect_equal(H2, test.hess21)


