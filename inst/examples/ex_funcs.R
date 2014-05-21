## ex_funcs.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun

## Functions to compute objective function, gradient, and Hessian for
## example.R


log.f <- function(pars, Y, X, inv.Omega, inv.Sigma, ...) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=FALSE)
  mu <- pars[(N*k+1):(N*k+k)]

  bx <- colSums(X * beta)
  
  log.p <- bx - log1p(exp(bx))
  log.p1 <- -log1p(exp(bx))
  
  data.LL <- sum(Y*log.p + (T-Y)*log.p1)
  
  Bmu <- apply(beta, 2, "-", mu)
  
  prior <- -0.5 * sum(diag(tcrossprod(Bmu) %*% inv.Sigma))
  hyp <- -0.5 * t(mu) %*% inv.Omega %*% mu
  res <- data.LL + prior + hyp
  return(as.numeric(res))
  
}

dlog.f.db <- function(pars, Y, X, inv.Omega, inv.Sigma) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=FALSE)
  mu <- pars[(N*k+1):length(pars)]
  bx <- colSums(X * beta)

  p <- exp(bx)/(1+exp(bx))

  tmp <- Y - T*p

  dLL.db <- apply(X,1,"*",tmp)
    
  Bmu <- apply(beta, 2, "-", mu)
  dprior <- -inv.Sigma %*% Bmu
  
  res <- t(dLL.db) + dprior

  return(as.vector(res))
 
}

dlog.f.dmu <- function(p, Y, X, inv.Omega, inv.Sigma) {

  beta <- matrix(p[1:(N*k)], k, N, byrow=FALSE)
  mu <- p[(N*k+1):length(p)]
  Bmu <- apply(beta, 2, "-", mu)

  res <- inv.Sigma %*% (rowSums(Bmu)) -  inv.Omega %*% mu
  return(res)
}

get.grad <- function(p, Y, X, inv.Omega, inv.Sigma, ...) {

  q1 <- dlog.f.db(p, Y, X, inv.Omega, inv.Sigma)
  q2 <- dlog.f.dmu(p, Y, X, inv.Omega, inv.Sigma)
  res <- c(q1, q2)
  return(res)
}


d2.db <- function(pars, Y, X, inv.Sigma, ...) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=FALSE)
  mu <- pars[(N*k+1):length(pars)]
  ebx <- exp(colSums(X * beta))

  p <- ebx/(1+ebx)
  
  q <- vector("list",length=N)
  for (i in 1:N) {
    q[[i]] <- -T*p[i]*(1-p[i])*tcrossprod(X[,i]) - inv.Sigma
##    q[[i]] <- -T*p[i]*(1-p[i])*XX.list[[i]] - inv.Sigma
  }
  B <- bdiag(q)
  return(B)
  
}

d2.dmu <- function(inv.Sigma, inv.Omega) {
  return(-N*inv.Sigma-inv.Omega)
}

d2.cross <- function(inv.Sigma) {
  res <- kronecker(Matrix(rep(1,N),nrow=1),inv.Sigma)
  return(res)
}


get.hess <- function(p, Y, X, inv.Omega, inv.Sigma, ...) {
  SX <- Matrix(inv.Sigma)
  XO <- Matrix(inv.Omega)
  B1 <- d2.db(p, Y, X, SX, ...)
  cross <- d2.cross(SX)
  Bmu <- d2.dmu(SX, XO)
  res <- rBind(cBind(B1,t(cross)),cBind(cross, Bmu))

  return(res)
}



get.hess.struct <- function(N, k) {

  B1 <- kronecker(Diagonal(N),Matrix(TRUE,k,k))
  B2 <- Matrix(TRUE,k,N*k)
  B3 <- Matrix(TRUE,k,k)
  H <- cBind(rBind(B1,B2),rBind(t(B2),B3))
  res <- Matrix.to.Coord(H)
  return(res)
  
}

