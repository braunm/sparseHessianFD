## example.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun.  See LICENSE file for details.

## this file demonstrates how to create a sparseHessianObj object, and use
## it to compute the Hessian of the function, given the gradient and sparsity
## structure.

print("Cleanup")
rm(list=ls())
print("gc")
gc()
print("moving on")

library(sparseHessianFD)
library(mvtnorm)
library(plyr)
library(Rcpp)
library(RcppEigen)

##source("ex_funcs.R")
print("set.seed")
set.seed(123)

N <- 5
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


hess.struct <- get.hess.struct(N, k)  ## for SparseFD method only

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


print("new obj")
obj <- new("sparseHessianFD", nvars, get.f, get.df) 
print("hess init")
obj$hessian.init(hess.struct$iRow, hess.struct$jCol, 0, 1e-5)
print("gc1")
gc()

## obj <- new.sparse.hessian.obj(par, log.f, get.grad, hess.struct, 
##                               Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)

print("Computing fn")
fn <- obj$fn(par)
print("gc2")
gc()
print("Computing gr")
gr <- obj$gr(par)
print("gc3")
gc()
print("Computing Hessian")
hs <- obj$hessian(par)
stop()
print("gc4")
gc()
print("Computing fn and gr")
fdf <- obj$fngr(par)
print("Done with obj")


print("rm(obj)")
rm(obj)
#gc()
print("rm(get.f)")
rm(get.f)
gc()
print("rm(get.df)")
rm(get.df)
gc()
print("end")

## print("Computing fn")
## fn1 <- get.fn(par, obj)
## print("Computing gr")
## gr1 <- get.gr(par, obj)
## print("Computing Hessian, again")
## hs1 <- get.hessian(par, obj)
## print("Computing fn and gr, again")
## get.fngr(par, obj)
## print("Last thing")
## ##H2 <- get.hess(par,  Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
## print("Done")



 

