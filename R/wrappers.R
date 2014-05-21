## wrappers.R --   Part of the sparsehessianFD package for the R programming language.
##
## Copyright (C) 2013 Michael Braun




new.sparse.hessian.obj <- function(x, fn, gr, hs, fd.method=0L, eps=sqrt(.Machine$double.eps), ...) {

  nvars <- length(x)
  
  if (!is.function(fn)) stop ("Error in new.sparse.hessian.obj:  fn must be a function")
  if (!is.function(gr)) stop ("Error in new.sparse.hessian.obj:  gr must be a function")
  if (is.null(hs)) stop("Error in new.sparse.hessian.obj: you must supply structure of the Hessian.")
  if (!is.list(hs)) stop ("Error in new.sparse.hessian.obj:  hs must be a list")
  if (!all.equal(names(hs),c("iRow","jCol"))) stop ("Error in new.sparse.hessian.obj.  Names of hs must be iRow and jCol")
  
  if (!is.integer(fd.method) || fd.method<0 || fd.method>1 || !is.finite(fd.method)) {
    stop("Error in get.sparse.hessian:  fd.method must be an integer, and either 0 or 1.")
  }

  fn1 <- function(x) fn(x,...)  ## currying the data in
  gr1 <- if (!is.null(gr)) function(x) gr(x,...)

  ## test fn and gr

  r1 <- fn1(x)
  if (!is.finite(r1)) stop("Error in new.sparse.hessian.obj:  fn at starting values is not finite.")
  r2 <- gr1(x)
  if (any(!is.finite(r2))) stop("Error in new.sparse.hessian.obj:  at least one element of gr is not finite.")
  if (length(r2)!=length(x)) stop("Error in new.sparse.hessian.obj:  length of gradient does not match length of starting values.")

  obj <- new("sparseHessianObj", nvars, fn1, gr1, hs, fd.method, eps)

  return(obj)

}

get.fn <- function(x, obj) obj$fn(x)

get.gr <- function(x, obj) obj$gr(x)

get.hessian <- function(x, obj) obj$hessian(x)

get.fngr <- function(x, obj) obj$fngr(x)
