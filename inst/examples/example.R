## example.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun.  See LICENSE file for details.

## this file demonstrates how to create a sparseHessianObj object, and use
## it to compute the Hessian of the function, given the gradient and sparsity
## structure.


library(plyr)
library(Rcpp)
library(RcppEigen)
library(sparseHessianFD)

source("ex_funcs.R")

set.seed(123)

N <- 15
k <- 3
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


hess.struct <- get.hess.struct(N, k)  ## for SparseFD method only

obj <- new.sparse.hessian.obj(par, log.f, get.grad, hess.struct, 
                              Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)

fn <- obj$fn(par)
gr <- obj$gr(par)
hs <- obj$hessian(par)
fdf <- obj$fngr(par)
fdfh <- obj$all(par)

fn1 <- get.fn(par, obj)
gr1 <- get.gr(par, obj)
hs1 <- get.hessian(par, obj)
get.fngr(par, obj)

H2 <- get.hess(par,  Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)




 

