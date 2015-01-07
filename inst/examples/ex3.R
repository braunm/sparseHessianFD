## example.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun.  See LICENSE file for details.

## this file demonstrates how to create a sparseHessianObj object, and use
## it to compute the Hessian of the function, given the gradient and sparsity
## structure.

rm(list=ls())
gc()


library(mvtnorm)
library(Rcpp)
library(RcppEigen)
library(numDeriv)



set.seed(123)


data(binary)
Y <- binary$Y
X <- binary$X
T <- binary$T
N <- length(Y)
k <- NROW(X)

nvars <- as.integer(N*k + k)
par <- rnorm(nvars) ## random starting values

Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)


make.f <- function(Y, X, inv.Omega, inv.Sigma) {
    function(pars) {
        log.f(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
    }
}

make.df <- function(Y, X, inv.Omega, inv.Sigma) {
    function(pars) {
        get.grad(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
    }
}

make.hess <- function(Y, X, inv.Omega, inv.Sigma) {
    function(pars) {
        get.hess(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
    }
}

get.f <- make.f(Y, X, inv.Omega, inv.Sigma)
get.df <- make.df(Y, X, inv.Omega, inv.Sigma)
get.d2f <- make.hess(Y, X, inv.Omega, inv.Sigma)

## True values for test

true.f <- get.f(par)
true.df <- get.df(par)
true.hess <- get.d2f(par)


## Get hessian structure
hess.struct <- Matrix.to.Coord(as(tril(true.hess), "lMatrix"))

obj <- new("sparseHessianFD", nvars, get.f, get.df) )
obj$hessian.init(hess.struct$iRow, hess.struct$jCol, 1, 1e-6)

## obj <- new.sparse.hessian.obj(par, log.f, get.grad, hess.struct, 
##                               Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)

fn <- obj$fn(par)
gr <- obj$gr(par)
hs <- obj$hessian(par)



 

