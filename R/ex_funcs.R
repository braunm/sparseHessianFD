## ex_funcs.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2013 Michael Braun

## Functions to compute objective function, gradient, and Hessian for
## binary choice example, and some unit tests.


#' @name binary
#' @title Binary choice example
#' @description Functions for binary choice example in the vignette.
#' @param P Numeric vector of length (N+1)*k.  First N*k elements are heterogeneous coefficients. The remaining k elements are population parameters.
#' @param Y,X Data matrices
#' @param inv.Omega,inv.Sigma Prior parameters
#' @return Log posterior density, gradient and Hessian.
#' @details Hessian is sparse, and returned as a dgcMatrix object
#' @rdname binary
binary.f <- function(P, Y, X, T, inv.Omega, inv.Sigma, order.row=FALSE) {

    N <- length(Y)
    k <- NROW(X)
    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):(N*k+k)]
    
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

#' @rdname binary
binary.grad <- function(P, Y, X, T, inv.Omega, inv.Sigma, order.row=FALSE) {

  q1 <- .dlog.f.db(P, Y, X, T, inv.Omega, inv.Sigma, order.row=order.row)
  q2 <- .dlog.f.dmu(P, Y, X, T, inv.Omega, inv.Sigma, order.row=order.row)
  res <- c(q1, q2)
  return(res)
}

#' @rdname binary
binary.hess <- function(P, Y, X, T, inv.Omega, inv.Sigma,  order.row=FALSE) {
    N <- length(Y)
    k <- NROW(X)
    
    
    SX <- Matrix(inv.Sigma)
    XO <- Matrix(inv.Omega)
    B1 <- .d2.db(P, Y, X, T, SX, order.row)
    
    if(order.row) {
        pvec <- NULL
        for (i in 1:N) {
            pvec <- c(pvec, seq(i, i+(k-1)*N, length=k))
        }
        
        pmat <- as(pvec, "pMatrix")
        B2 <- Matrix::t(pmat) %*% B1 %*% pmat
        cross <- .d2.cross(N, SX) %*% pmat      
    } else {
        B2 <- B1
        cross <- .d2.cross(N, SX)
    }
    
    Bmu <- .d2.dmu(N,SX, XO)
    res <- rBind(cBind(B2, Matrix::t(cross)),cBind(cross, Bmu))
    
    return(res)
}



.dlog.f.db <- function(P, Y, X, T, inv.Omega, inv.Sigma, order.row) {

    N <- length(Y)
    k <- NROW(X)
    
    
    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
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

.dlog.f.dmu <- function(P, Y, X, T, inv.Omega, inv.Sigma, order.row) {

    N <- length(Y)
    k <- NROW(X)
    
    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
    Bmu <- apply(beta, 2, "-", mu)
    
    res <- inv.Sigma %*% (rowSums(Bmu)) -  inv.Omega %*% mu
    return(res)
}



.d2.db <- function(P, Y, X, T, inv.Sigma, order.row) {

    N <- length(Y)
    k <- NROW(X)
    
    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
    ebx <- exp(colSums(X * beta))
    
    p <- ebx/(1+ebx)
    
    q <- vector("list",length=N)
    for (i in 1:N) {
        q[[i]] <- -T*p[i]*(1-p[i])*tcrossprod(X[,i]) - inv.Sigma
    }
    B <- bdiag(q)
    return(B)    
}


.d2.dmu <- function(N, inv.Sigma, inv.Omega) {
  return(-N*inv.Sigma-inv.Omega)
}

.d2.cross <- function(N, inv.Sigma) {
  res <- kronecker(Matrix(rep(1,N),nrow=1),inv.Sigma)
  return(res)
}




