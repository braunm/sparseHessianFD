## ----setup1, echo=FALSE, cache=FALSE-------------------------------------
suppressPackageStartupMessages(require(Matrix))

## ----setup2, echo=FALSE--------------------------------------------------
N <- 5
k <- 2
nv1 <- (N + 1) * k
nels1 <- nv1^2
nnz1 <- (N + 1) * k^2 + 2 * N * k^2
nnz1LT <- (N + 1) * k * (k+1)/2 + N * k^2
Q <- 1000
nv2 <- (Q+1) * k
nels2 <- nv2^2
nnz2 <- (Q+1) * k^2 + 2 * Q * k^2
nnz2LT <- (Q+1) * k * (k+1)/2 + Q * k^2


## ----blockarrow, echo=FALSE----------------------------------------------
Mat <- as(kronecker(diag(N), matrix(1, k, k)),"sparseMatrix")
Mat <- rBind(Mat, Matrix(1, k, N*k))
Mat <- cBind(Mat, Matrix(1, k * (N+1), k))
cat("\nFigure 1a:  block-arrow pattern\n")
printSpMatrix(as(Mat,"nMatrix"))

## ----banded, echo=FALSE--------------------------------------------------
Mat <- kronecker(Matrix(1, k, k), diag(N))
Mat <- rBind(Mat, Matrix(1, k, N * k))
Mat <- cBind(Mat, Matrix(1, k * (N + 1), k))
cat("\nFigure 1b:  banded pattern\n")
printSpMatrix(as(Mat,"nMatrix"))

## ------------------------------------------------------------------------
library("sparseHessianFD")
bd <- kronecker(diag(3), matrix(TRUE, 2, 2))
Mat <- as(bd, "nMatrix")
cat("\nPatterns from page 10\n")
printSpMatrix(tril(Mat))
mc <- Matrix.to.Coord(tril(Mat))
mc

## ------------------------------------------------------------------------
pattern <- sparseMatrix(i=mc$rows, j=mc$cols)
printSpMatrix(pattern)


## ----binaryInit, tidy=FALSE----------------------------------------------
set.seed(123)
data("binary")
cat("\nstructure of data frame with sample data\n")
str(binary)
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
T <- binary[["T"]]
nvars <- as.integer(N * k + k)
priors <- list(inv.Sigma = rWishart(1, k + 5, diag(k))[,,1],
   inv.Omega = diag(k))

## ----trueValues----------------------------------------------------------
P <- rnorm(nvars)
order.row <- FALSE
true.f <- binary.f(P, binary, priors, order.row=order.row)
true.grad <- binary.grad(P, binary, priors, order.row=order.row)
true.hess <- binary.hess(P, binary, priors, order.row=order.row)

## ----binaryRowsCols------------------------------------------------------
pattern <- Matrix.to.Coord(tril(true.hess))
cat("\nstructure of output from Matrix.to.Coord function\n")
str(pattern)

## ----usingSparseHessianFD1, tidy=FALSE-----------------------------------
obj <- sparseHessianFD(P, fn=binary.f, gr=binary.grad,
   rows=pattern[["rows"]], cols=pattern[["cols"]],
   data=binary, priors=priors, order.row=order.row)
f <- obj$fn(P)
identical(f, true.f)
gr <- obj$gr(P)
identical(gr, true.grad)


## ----usingSparseHessianFD2-----------------------------------------------
hs <- obj$hessian(P)
cat("\nRelative error using finite difference method (p. 13)\n")
print(mean(abs(hs-true.hess))/mean(abs(hs)))

## ----tidy=FALSE----------------------------------------------------------
obj2 <- sparseHessianFD(P, fn=binary.f, gr=binary.grad,
   rows=pattern[["rows"]], cols=pattern[["cols"]],
   complex=TRUE,
   data=binary, priors=priors, order.row=order.row)
hs2 <- obj2$hessian(P)
cat("\nRelative error using complex step method (p. 13-14)\n")
print(mean(abs(hs2-true.hess))/mean(abs(hs2)))

