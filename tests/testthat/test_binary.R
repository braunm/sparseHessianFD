context("binary example")

test_that("binary_example", {
    require(numDeriv)

    set.seed(123)


    data(binary)
    Y <- binary$Y
    X <- binary$X
    T <- binary$T
    N <- length(Y)
    k <- NROW(X)
    nvars <- as.integer(N*k + k)
    P <- rnorm(nvars) ## random starting values

    Omega <- diag(k)
    inv.Omega <- solve(Omega)
    inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]

    make.funcs <- function(Y, X, inv.Omega, inv.Sigma) {
        res <- vector("list", length=3)
        names(res) <- c("fn", "gr", "hessian")
        res$fn <-  function(pars) {
            log.f(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
        }
        res$gr <-  function(pars) {
            get.grad(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
        }
        res$hessian <-  function(pars) {
            get.hess(pars, Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
        }
        return(res)
    }
        
    f <- make.funcs(Y, X, inv.Omega, inv.Sigma)
    
    ## True values for test
    
    true.val <- f$fn(P)
    true.grad <- f$gr(P)
    true.hess <- f$hessian(P)    
    
    ## Get hessian structure
    hess.struct <- Matrix.to.Coord(as(tril(true.hess), "lMatrix"))
    
    obj <- new("sparseHessianFD", nvars, f$fn, f$gr)
    obj$hessian.init(hess.struct$iRow, hess.struct$jCol, 1, 1e-6)
    
    ## obj <- new.sparse.hessian.obj(par, log.f, get.grad, hess.struct, 
    ##                               Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma)
    
    test.val <- obj$fn(P)
    test.grad <- obj$gr(P)
    test.hess <- obj$hessian(P)

    expect_equal(test.val, true.val)
    expect_equal(test.grad, true.grad)
    expect_equal(test.hess, true.hess)
})
    


 

