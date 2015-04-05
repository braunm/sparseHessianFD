library(Matrix)
library(testthat)

set.seed(123)
data(binary)
binary <- binary

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


time <- array(NA, dim=c(2,2,5,2))
dimnames(time) <- list(order=c("block","band"),
                       method=c("direct","indirect"),
                       stat=c("init","diff", "subst","fd","setup"),
                       version=c("old","new"))

time["block","indirect","init","new"] <- system.time(
    W1 <- color.cols(pattern1$rows, pattern1$cols))[["elapsed"]]
time["band","indirect","init","new"] <- system.time(
    W2 <- color.cols(pattern2$rows, pattern2$cols))[["elapsed"]]

delta <- 1e-7

time["block","indirect","diff","new"] <- system.time(
    Y1 <- get.diffs(P, df=f1$gr, pattern1$rows, pattern1$cols, W1, delta))[["elapsed"]]
time["band","indirect","diff","new"] <- system.time(
    Y2 <- get.diffs(P, df=f2$gr, pattern2$rows, pattern2$cols, W2, delta))[["elapsed"]]

time["block","indirect","subst","new"] <- system.time(
    H1 <- subst(Y1, W1, pattern1$rows, pattern1$cols, delta))[["elapsed"]]
time["band","indirect","subst","new"] <- system.time(
    H2 <- subst(Y2, W2, pattern2$rows, pattern2$cols, delta))[["elapsed"]]

time["block","indirect","setup","old"] <- system.time(
    obj10 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))[["elapsed"]]
time["block","indirect","init","old"] <- system.time(
    obj10$hessian.init(pattern1$rows, pattern1$cols, 0, delta))[["elapsed"]]

time["block","direct","setup","old"] <- system.time(
    obj11 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))[["elapsed"]]
time["block","direct","init","old"] <- system.time(
    obj11$hessian.init(pattern1$rows, pattern1$cols, 1, delta))[["elapsed"]]

time["band","indirect","setup","old"] <- system.time(
    obj20 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))[["elapsed"]]
time["band","indirect","init","old"] <- system.time(
    obj20$hessian.init(pattern2$rows, pattern2$cols, 0, delta))[["elapsed"]]

time["band","direct","setup","old"] <- system.time(
    obj21 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))[["elapsed"]]
time["band","direct","init","old"] <- system.time(
    obj21$hessian.init(pattern2$rows, pattern2$cols, 1, delta))[["elapsed"]]


time["block","indirect","fd","old"] <- system.time(test.hess10 <- obj10$hessian(P))[["elapsed"]]
time["block","direct","fd","old"] <- system.time(test.hess11 <- obj11$hessian(P))[["elapsed"]]
time["band","indirect","fd","old"] <- system.time(test.hess20 <- obj20$hessian(P))[["elapsed"]]
time["band","direct","fd","old"] <- system.time(test.hess21 <- obj21$hessian(P))[["elapsed"]]


expect_equal(true.hess1, H1, tolerance=1e-6)
expect_equal(true.hess2, H2, tolerance=1e-6)
expect_equal(H1, test.hess10)
expect_equal(H1, test.hess11)
expect_equal(H2, test.hess20)
expect_equal(H2, test.hess21)


