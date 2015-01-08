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

N <- 4
k <- 2
T <- 10
order.row <- TRUE

## Simulate data and set priors

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k, 0, 1)
Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)

X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- sapply(log.p, function(q) return(rbinom(1,T,exp(q))))

nvars <- as.integer(N*k + k)
par <- rnorm(nvars) ## random starting values
print("Creating true hessian")
true.hess <- drop0(binary.hess(par, Y, X, T, inv.Omega, inv.Sigma, order.row=order.row))
hess.struct <- Matrix.to.Coord(as(tril(true.hess), "lMatrix"))

make.f <- function(Y, X, T, inv.Omega, inv.Sigma, order.row=FALSE) {
    function(pars) {
        binary.f(pars, Y=Y, X=X, T=T, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma, order.row=order.row)
    }
}

make.df <- function(Y, X, T, inv.Omega, inv.Sigma, order.row=FALSE) {
    function(pars) {
        binary.grad(pars, Y=Y, X=X, T=T, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma, order.row=order.row)
    }
}

get.f <- make.f(Y, X, T, inv.Omega, inv.Sigma, order.row)
get.df <- make.df(Y, X, T, inv.Omega, inv.Sigma, order.row)


print("Creating object")
new.time <- system.time(obj <- new("sparseHessianFD", nvars, get.f, get.df) )
print(new.time)
print("Initializing")
init.time <- system.time(obj$hessian.init(hess.struct$iRow, hess.struct$jCol, 1, 1e-8))
print(init.time)


cat("Function eval:\n")
fn.time <- system.time(fn <- obj$fn(par))
print(fn.time)
cat("\nGradient eval:\n")
grad.time <- system.time(gr <- obj$gr(par))
print(grad.time)
cat("\nHessian eval:\n")
hess.time <- system.time(hs <- obj$hessian(par))
print(hess.time)
h2 <- drop0(hessian(get.f, par), 1e-8)
print("Correct result:?")

num.grad <- grad(get.f, par)
print(all.equal(gr, num.grad))
print(all.equal(hs, true.hess))
print(all.equal(hs, h2))




 

