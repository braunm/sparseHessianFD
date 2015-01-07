## example.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun.  See LICENSE file for details.

## this file demonstrates how to create a sparseHessianObj object, and use
## it to compute the Hessian of the function, given the gradient and sparsity
## structure.

rm(list=ls())
gc()
library(sparseHessianFD)
library(mvtnorm)
library(plyr)
library(Rcpp)
library(RcppEigen)
library(numDeriv)

set.seed(123)

N <- 3
k <- 2
T <- 20


## Simulate data and set priors

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k,0,10)
Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)
X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- laply(log.p, function(q) return(rbinom(1,T,exp(q))))

nvars <- as.integer(N*k + k)
par <- rnorm(nvars) ## random starting values


#hess.struct <- get.hess.struct(N, k)  ## for SparseFD method only

true.hess <- get.hess(par, Y, X, inv.Omega, inv.Sigma)
hess.struct <- Matrix.to.Coord(as(tril(true.hess), "lMatrix"))

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

get.f2 <- make.f(Y, X, inv.Omega, inv.Sigma)
get.df2 <- make.df(Y, X, inv.Omega, inv.Sigma)

get.f <- function(v) log.f(v, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
get.df <- function(v) get.grad(v, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)


obj <- new("sparseHessianFD", nvars, get.f2, get.df2) 
obj$hessian.init(hess.struct$iRow, hess.struct$jCol, 0, 1e-4)

## obj <- new.sparse.hessian.obj(par, log.f, get.grad, hess.struct, 
##                               Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)

fn <- obj$fn(par)
gr <- obj$gr(par)
hs <- obj$hessian(par)
h2 <- drop0(hessian(get.f, par), 1e-8)
fdf <- obj$fngr(par)

ha <- get.hess(par, Y, X, inv.Omega, inv.Sigma)




 

