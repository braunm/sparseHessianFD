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

    make.funcs <- function(Y, X, T, inv.Omega, inv.Sigma, order.row) {
        res <- vector("list", length=3)
        names(res) <- c("fn", "gr", "hessian")
        res$fn <-  function(pars) {
            binary.f(pars, Y=Y, X=X, T=T, inv.Omega=inv.Omega,
                     inv.Sigma=inv.Sigma, order.row=order.row)
        }
        res$gr <-  function(pars) {
            binary.grad(pars, Y=Y, X=X, T=T, inv.Omega=inv.Omega,
                        inv.Sigma=inv.Sigma, order.row=order.row)
        }
        res$hessian <-  function(pars) {
            binary.hess(pars, Y=Y, X=X, T=T, inv.Omega=inv.Omega,
                        inv.Sigma=inv.Sigma, order.row=order.row)
        }
        return(res)
    }
        
    f1 <- make.funcs(Y, X, T, inv.Omega, inv.Sigma, order.row=FALSE) ## block-arrow
    f2 <- make.funcs(Y, X, T, inv.Omega, inv.Sigma, order.row=TRUE) ## off-diagonals
    
    ## True values for test
    
    true.val1 <- f1$fn(P)
    true.grad1 <- f1$gr(P)
    true.hess1 <- f1$hessian(P)    

    true.val2 <- f2$fn(P)
    true.grad2 <- f2$gr(P)
    true.hess2 <- f2$hessian(P)    

    ## Get hessian structure
    hess.struct1 <- Matrix.to.Coord(tril(true.hess1))
    hess.struct2 <- Matrix.to.Coord(tril(true.hess2))
    
    obj10 <- new("sparseHessianFD", nvars, f1$fn, f1$gr)
    obj10$hessian.init(hess.struct1$iRow, hess.struct1$jCol, 0, 1e-7)

    obj11 <- new("sparseHessianFD", nvars, f1$fn, f1$gr)
    obj11$hessian.init(hess.struct1$iRow, hess.struct1$jCol, 1, 1e-7)

    obj20 <- new("sparseHessianFD", nvars, f2$fn, f2$gr)
    obj20$hessian.init(hess.struct2$iRow, hess.struct2$jCol, 0, 1e-7)

    obj21 <- new("sparseHessianFD", nvars, f2$fn, f2$gr)
    obj21$hessian.init(hess.struct2$iRow, hess.struct2$jCol, 1, 1e-7)


    test.val10 <- obj10$fn(P)
    test.grad10 <- obj10$gr(P)
    test.hess10 <- obj10$hessian(P)
    
    test.val11 <- obj11$fn(P)
    test.grad11 <- obj11$gr(P)
    test.hess11 <- obj11$hessian(P)
    
    test.val20 <- obj20$fn(P)
    test.grad20 <- obj20$gr(P)
    test.hess20 <- obj20$hessian(P)
    
    test.val21 <- obj21$fn(P)
    test.grad21 <- obj21$gr(P)
    test.hess21 <- obj21$hessian(P)

    expect_equal(test.val10, true.val1)
    expect_equal(test.grad10, true.grad1)
    expect_equal(test.hess10, true.hess1)
    expect_equal(test.val11, true.val1)
    expect_equal(test.grad11, true.grad1)
    expect_equal(test.hess11, true.hess1)
    expect_equal(test.val20, true.val2)
    expect_equal(test.grad20, true.grad2)
    expect_equal(test.hess20, true.hess2)
    expect_equal(test.val21, true.val2)
    expect_equal(test.grad21, true.grad2)
    expect_equal(test.hess21, true.hess2)
})
    


 

