library(sparseHessianFD)
library(mvtnorm)
library(plyr)
library(microbenchmark)
library(doParallel)
library(numDeriv)

run.par <- TRUE
if (run.par) registerDoParallel(cores=12) else registerDoParallel(cores=1)

binary_sim <- function(N, k, T) {
    x.mean <- rep(0,k)
    x.cov <- diag(k)
    x.cov[1,1] <- .02
    x.cov[k,k] <- x.cov[1,1]
    mu <- seq(-2,2,length=k)
    Omega <- diag(k)

    X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
    B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
    XB <- colSums(X * B)
    log.p <- XB - log1p(exp(XB))
    Y <- sapply(log.p, function(q) return(rbinom(1,T,exp(q))))
    list(Y=Y, X=X, T=T)
}

priors_sim <- function(k) {
    list(inv.Omega = solve(diag(k)),
         inv.Sigma = rWishart(1,k+5,diag(k))[,,1])
}


make_funcs <- function(D, priors) {
    res <- vector("list", length=3)
    names(res) <- c("fn", "gr", "hessian")
    res$fn <-  function(pars) binary.f(pars, data=D, priors=priors)
    res$gr <-  function(pars) binary.grad(pars, data=D, priors=priors)
    res$hessian <-  function(pars) binary.hess(pars, data=D, priors=priors)
    return(res)
}

run_test <- function(Nk, reps=50) {
    N <- as.numeric(Nk["N"])
    k <- as.numeric(Nk["k"])
    data <- binary_sim(N, k, T=20)
    priors <- priors_sim(k)
    F <- make_funcs(D=data, priors=priors)
    nvars <- N*k+k
    M <- as(Matrix::kronecker(Matrix::Diagonal(N),Matrix(1,k,k)),"nMatrix") %>%
      rBind(Matrix(TRUE,k,N*k)) %>%
      cBind(Matrix(TRUE, k*(N+1), k)) %>%
      as("nMatrix")
    pat <- Matrix.to.Coord(tril(M))
    X <- rnorm(nvars)
    obj <- sparseHessianFD(X, F$fn, F$gr, pat$rows, pat$cols)

    bench <- microbenchmark(
        numDeriv = hessian(obj$fn, X),
        df = obj$gr(X),
        sparse = obj$hessian(X))
    vals <- ddply(data.frame(bench), "expr",
                  function(x) return(data.frame(expr=x$expr,
                                                time=x$time,
                                                rep=1:length(x$expr))))
    res <- data.frame(N=N, k=k,
                      bench=vals)
    cat("Completed N = ",N,"\tk = " , k ,"\n")
    return(res)
}


cases <- expand.grid(k=c(2,3,4),
                     N=c(9, 12, 15))

runs <- adply(cases, 1, run_test, reps=20, .parallel=run.par)

tab4 <-  mutate(runs, ms=bench.time/1000000) %>%
  select(-bench.time) %>%
  spread(bench.expr, ms) %>%
  gather(method, hessian, c(numDeriv, sparse)) %>%
  mutate(M=N*k+k, hessian.df=hessian/df) %>%
  gather(stat, time, c(hessian, hessian.df)) %>%
  group_by(N, k, method, M, stat)  %>%
  summarize(mean=mean(time), sd=sd(time)) %>%
  gather(stat2, value, mean:sd) %>%
  dcast(N+k+M~stat+method+stat2,value.var="value") %>%
  arrange(M)







