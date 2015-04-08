rm(list=ls())
gc()
library(testthat)

set.seed(123)
data(binary_super)
binary <- binary_super

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

## Get hessian structure2
pattern1LT <- Matrix.to.Coord(tril(true.hess1))
pattern2LT <- Matrix.to.Coord(tril(true.hess2))
pattern1 <- Matrix.to.Coord(true.hess1)
pattern2 <- Matrix.to.Coord(true.hess2)

time <- array(NA, dim=c(2,2,3,2))
dimnames(time) <- list(order=c("block","band"),
                       method=c("direct","indirect"),
                       stat=c("init","fd","setup"),
                       version=c("old","new"))


rows1 <- pattern1$rows
rows2 <- pattern2$rows
cols1 <- pattern1$cols
cols2 <- pattern2$cols

rows1LT <- pattern1LT$rows
rows2LT <- pattern2LT$rows
cols1LT <- pattern1LT$cols
cols2LT <- pattern2LT$cols

time["block","indirect","init","new"] <- system.time(
    W1 <- color.cols(rows1LT, cols1LT))[["elapsed"]]
time["band","indirect","init","new"] <- system.time(
    W2 <- color.cols(rows2LT, cols2LT))[["elapsed"]]

delta <- 1e-7

time["block","indirect","fd","new"] <- system.time(
    H1 <- get.fd(P, df=f1$gr, rows1LT, cols1LT, W1, delta))[["elapsed"]]
time["band","indirect","fd","new"] <- system.time(
    H2 <- get.fd(P, df=f2$gr, rows2LT, cols2LT, W2, delta))[["elapsed"]]

time["block","indirect","setup","old"] <- system.time(
    obj10 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))[["elapsed"]]
time["block","indirect","init","old"] <- system.time(
    obj10$hessian.init(rows1LT, cols1LT, 0, delta))[["elapsed"]]

time["block","direct","setup","old"] <- system.time(
    obj11 <- new("sparseHessianFD", nvars, f1$fn, f1$gr))[["elapsed"]]
time["block","direct","init","old"] <- system.time(
    obj11$hessian.init(rows1LT, cols1LT, 1, delta))[["elapsed"]]

time["band","indirect","setup","old"] <- system.time(
    obj20 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))[["elapsed"]]
time["band","indirect","init","old"] <- system.time(
    obj20$hessian.init(rows2LT, cols2LT, 0, delta))[["elapsed"]]

time["band","direct","setup","old"] <- system.time(
    obj21 <- new("sparseHessianFD", nvars, f2$fn, f2$gr))[["elapsed"]]
time["band","direct","init","old"] <- system.time(
    obj21$hessian.init(rows2LT, cols2LT, 1, delta))[["elapsed"]]


time["block","indirect","fd","old"] <- system.time(test.hess10 <- obj10$hessian(P))[["elapsed"]]
time["block","direct","fd","old"] <- system.time(test.hess11 <- obj11$hessian(P))[["elapsed"]]
time["band","indirect","fd","old"] <- system.time(test.hess20 <- obj20$hessian(P))[["elapsed"]]
time["band","direct","fd","old"] <- system.time(test.hess21 <- obj21$hessian(P))[["elapsed"]]


expect_equivalent(true.hess1, drop0(H1))
expect_equivalent(true.hess2, drop0(H2))
expect_equivalent(drop0(H1), test.hess10)
expect_equivalent(drop0(H1), test.hess11)
expect_equivalent(drop0(H2), test.hess20)
expect_equivalent(drop0(H2), test.hess21)


