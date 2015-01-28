## ----, echo=FALSE--------------------------------------------------------
require(Matrix)
require(sparseHessianFD)
N <- 6
k <- 2
nv1 <- (N+1)*k
nels1 <- nv1^2
nnz1 <- (N+1)*k^2 + 2*N*k^2
nnz1LT <- (N+1)*k*(k+1)/2 + N*k^2
Q <- 1000
nv2 <- (Q+1)*k
nels2 <- nv2^2
nnz2 <- (Q+1)*k^2 + 2*Q*k^2
nnz2LT <- (Q+1)*k*(k+1)/2 + Q*k^2
options(scipen=999)

## ----, echo=FALSE--------------------------------------------------------
M <- as(kronecker(diag(N),matrix(1,k,k)),"lMatrix")
M <- rBind(M, Matrix(TRUE,k,N*k))
M <- cBind(M, Matrix(TRUE, k*(N+1), k))
print(M)

## ----, echo=FALSE--------------------------------------------------------
M <- as(kronecker(matrix(1,k,k), diag(N)),"lMatrix")
M <- rBind(M, Matrix(TRUE,k,N*k))
M <- cBind(M, Matrix(TRUE, k*(N+1), k))
print(M)

## ----, collapse=TRUE-----------------------------------------------------
set.seed(123)
data(binary)
str(binary)
N <- length(binary$Y)
k <- NROW(binary$X)
nvars <- as.integer(N*k + k)
P <- rnorm(nvars) ## random starting values
priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))


## ----, echo=FALSE--------------------------------------------------------
options(scipen=-999)

## ----, tidy=TRUE---------------------------------------------------------
true.f <- binary.f(P, binary, priors, order.row=FALSE)
true.grad <- binary.grad(P, binary, priors, order.row=FALSE)
true.hess <- binary.hess(P, binary, priors, order.row=FALSE)

## ----, collapse=TRUE-----------------------------------------------------
pattern <- Matrix.to.Coord(tril(true.hess))
str(pattern)	    	    

## ------------------------------------------------------------------------
obj <- sparseHessianFD.new(P, binary.f, binary.grad,
       rows=pattern$rows, cols=pattern$cols,
       data=binary, priors=priors,
       order.row=FALSE)

## ------------------------------------------------------------------------
f <- obj$fn(P)
gr <- obj$gr(P)
hs <- obj$hessian(P)

## ----, collapse=TRUE-----------------------------------------------------


all.equal(f, true.f)
all.equal(gr, true.grad)
all.equal(hs, true.hess)	      	     	     

## ----, echo=FALSE--------------------------------------------------------
options(scipen=0)

## ----, tidy=TRUE---------------------------------------------------------
library(numDeriv)
hess.time <- system.time(H1 <- obj$hessian(P))	  
fd.time <- system.time(H2 <- hessian(obj$fn, P))
H2 <- drop0(H2, 1e-7) ## treat values < 1e-8 as zero
print(hess.time)
print(fd.time)	   

