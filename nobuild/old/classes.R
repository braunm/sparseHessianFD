## classes.R --   Part of the sparsehessianFD package for the R programming language.
##
## Copyright (C) 2013 Michael Braun


setClass("sparseHessianObj", representation(pointer = "externalptr") )

.sparseHessianObj_method <- function(name) {
  paste("sparseHessianObj",name,sep="__")
}


setMethod("$", "sparseHessianObj", function(x, name) {
  function(...) .Call(.sparseHessianObj_method(name), x@pointer, ...)
})


setMethod("initialize", "sparseHessianObj", function(.Object, nvars, fn, gr, hs, fd_method, eps) {
  .Object@pointer <- .Call(.sparseHessianObj_method("new"), nvars, fn, gr, hs, fd_method, eps)
  .Object
} )




