## ex_funcs.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun

## Functions to compute objective function, gradient, and Hessian for
## example.R


log.f <- function(pars, Y, X, inv.Omega, inv.Sigma, ..., order.row=FALSE) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=order.row)
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


get.grad <- function(p, ..., order.row=FALSE) {

  q1 <- dlog.f.db(p, ..., order.row=order.row)
  q2 <- dlog.f.dmu(p, ..., order.row=order.row)
  res <- c(q1, q2)
  return(res)
}


get.hess <- function(p, Y, X, inv.Omega, inv.Sigma, ...., order.row=FALSE) {
  SX <- Matrix(inv.Sigma)
  XO <- Matrix(inv.Omega)
  B1 <- d2.db(p, Y, X, SX, order.row)

  if(order.row) {
      pvec <- NULL
      for (i in 1:N) {
          pvec <- c(pvec, seq(i, i+(k-1)*N, length=k))
      }
      
      pmat <- as(pvec, "pMatrix")
      B2 <- Matrix::t(pmat) %*% B1 %*% pmat
      cross <- d2.cross(SX) %*% pmat      
  } else {
      B2 <- B1
      cross <- d2.cross(SX)
  }
  

  Bmu <- d2.dmu(SX, XO)
  res <- rBind(cBind(B2, Matrix::t(cross)),cBind(cross, Bmu))

  return(res)
}



dlog.f.db <- function(pars, Y, X, inv.Omega, inv.Sigma, order.row) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=order.row)
  mu <- pars[(N*k+1):length(pars)]
  bx <- colSums(X * beta)

  p <- exp(bx)/(1+exp(bx))

  tmp <- Y - T*p

  dLL.db <- apply(X,1,"*",tmp)
    
  Bmu <- apply(beta, 2, "-", mu)
  dprior <- -inv.Sigma %*% Bmu
  
  res <- t(dLL.db) + dprior
  if (order.row) res <- t(res)

  return(as.vector(res))
 
}

dlog.f.dmu <- function(p, Y, X, inv.Omega, inv.Sigma, order.row) {

  beta <- matrix(p[1:(N*k)], k, N, byrow=order.row)
  mu <- p[(N*k+1):length(p)]
  Bmu <- apply(beta, 2, "-", mu)

  res <- inv.Sigma %*% (rowSums(Bmu)) -  inv.Omega %*% mu
  return(res)
}



d2.db <- function(pars, Y, X, inv.Sigma, order.row) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=order.row)
  mu <- pars[(N*k+1):length(pars)]
  ebx <- exp(colSums(X * beta))

  p <- ebx/(1+ebx)
  
  q <- vector("list",length=N)
  for (i in 1:N) {
    q[[i]] <- -T*p[i]*(1-p[i])*tcrossprod(X[,i]) - inv.Sigma
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




