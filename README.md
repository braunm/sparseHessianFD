---
title:  The sparseHessianFD package
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# Why use the sparseHessianFD package?

Suppose you have a scalar-valued objective function of _N_ variables.  You have an **R** function that returns the function value, and another that returns the gradient.  You still need the
Hessian, which would be hard to derive analytically, and translate
into an **R** function.  If _N_ is very large, it can take a long time
to estimate the Hessian numerically, and a lot of memory to store all _N^2_
elements.

This scalability problem comes up in at least two areas of statistics:

1.  using derivative-based algorithms for nonlinear optimization; and
2.  sampling from, and estimating the density of, high-dimensional
multivariate normal distributions.

The *sparseHessianFD* package solves this problem when the Hessian is
sparse, and the sparsity pattern is known in advance.  The sparsity pattern is just the row and column indices of the
non-zero elements in the lower triangle of the Hessian.  Sparse Hessians arise, for example, in hierarchical models in
which parameters are conditionally independent across heterogeneous
units.  In that case, we know which cross-partial derivatives are zero.

See the vignette for more about sparsity patterns in Hessians.


# Getting the package


The package is available on CRAN, so a simple

```
install.packages("sparseHessianFD")
```

should suffice.  The package is being developed on github at
github.com/braunm/sparseHessian, and you can get development branches
there.

# Using the package

To use the package, you need at least the following:

-  an R function that takes a parameter vector as the first argument,
   and returns the value of the function for which you need the
   Hessian;
-  an R function that takes the same arguments as above, and returns
the gradient of the target function.
-  two integer vectors, one for the row indices of the non-zero
   elements of the lower triangle of the Hessian, and one for the
   column indices.

The `sparseHessianFD` function returns an object of class
sparseHessianFD.  This class contains methods to return the function
value, gradient and Hessian.  The Hessian is returned in a sparse,
compressed format as a `dgCMatrix` object (defined in the *Matrix*
package).

A typical usage might look like this.


```r
P <- rnorm(N) ## random "starting value" for initialization
obj <- sparseHessianFD(P, f, df, rows, cols)
X <- rnorm(N) ## another set of variables at which the function is evaluated
val <- obj$fn(X) ## function value
gr <- obj$gr(X) ## gradient
hess <- obj$hessian(X)  ## Hessian, as a sparse dgCMatrix object

```

See the vignette and package manual for more details, options,
convenience functions, etc.
