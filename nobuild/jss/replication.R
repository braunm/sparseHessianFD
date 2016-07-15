library("sparseHessianFD")
library("Matrix")
library("mvtnorm")
library("microbenchmark")
library("doParallel")
library("dplyr")
library("tidyr")
library("numDeriv")
library("reshape2")
library("ggplot2")


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
  res$fn <-  function(pars) binary.f(pars, data=D, priors=priors, order.row=TRUE)
  res$gr <-  function(pars) binary.grad(pars, data=D, priors=priors, order.row=TRUE)
  res$hessian <-  function(pars) binary.hess(pars, data=D, priors=priors, order.row=TRUE)
  return(res)
}

make_perm <- function(iRow, jCol) {
  tmp <- sparseMatrix(i=iRow, j=jCol, index1=TRUE, symmetric=TRUE)
  perm <- order(Matrix::rowSums(tmp), decreasing=TRUE)
  invperm <- invPerm(perm)
  L <- tril(tmp[perm,perm])
  ptr <- Matrix.to.Pointers(L, order="row", index1=TRUE)
  idx <- ptr[[1]]
  pntr <- ptr[[2]]
}

get_id <- function(x) 1:length(x)

run_test_fig4 <- function(NkT, reps=50) {
  ## Replication function for Figure 4
  N <- as.numeric(NkT["N"])
  k <- as.numeric(NkT["k"])
  T <- as.numeric(NkT["T"])
  data <- binary_sim(N, k, T)
  priors <- priors_sim(k)
  F <- make_funcs(D=data, priors=priors)
  nvars <- N*k+k

  M <- as(Matrix::kronecker(Matrix(1,k,k), Matrix::Diagonal(N)),"nMatrix") %>%
    rBind(Matrix(TRUE,k,N*k)) %>%
    cBind(Matrix(TRUE, k*(N+1), k)) %>%
    as("nMatrix")

  pat <- Matrix.to.Coord(tril(M))
  X <- rnorm(nvars)

  obj <- sparseHessianFD(X, F$fn, F$gr, pat$rows, pat$cols)
  colors <- obj$partition()
  perm <- obj$get_perm()
  ncolors <- length(unique(colors))
  nvars <- N*k+k

  bench <- microbenchmark(
      setup = sparseHessianFD(X, F$fn, F$gr, pat$rows, pat$cols),
      colors = coloring(tril(M[perm,perm])),
      perm = make_perm(pat$rows, pat$cols),
      f = obj$fn(X),
      df = obj$gr(X),
      hess = obj$hessian(X),
      times = reps)

  vals <- plyr::ddply(data.frame(bench), "expr",
                function(x) return(data.frame(expr=x$expr,
                                              time=x$time,
                                              rep=1:length(x$expr))))

  res <- data.frame(N=N, k=k, T=T,
                    bench=vals,
                    ncolors=ncolors)

  cat("Completed N = ",N,"\tk = " , k , "\tT = ",T,"\n")

  return(res)

}


run_test_tab4 <- function(Nk, reps=50) {
    ## Replication function for Table 4
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
        numDeriv = numDeriv::hessian(obj$fn, X),
        df = obj$gr(X),
        sparse = obj$hessian(X))
    vals <- plyr::ddply(data.frame(bench), "expr",
                  function(x) return(data.frame(expr=x$expr,
                                                time=x$time,
                                                rep=1:length(x$expr))))
    res <- data.frame(N=N, k=k,
                      bench=vals)
    cat("Completed N = ",N,"\tk = " , k ,"\n")
    return(res)
}

## Replicate Figure 4

cases_fig4 <- expand.grid(k=c(8, 6, 4, 2),
                          N=c(25, 50, seq(75,2500, by=75)),
                          T=20
                          )

runs_fig4 <- plyr::adply(cases_fig4, 1, run_test_fig4, reps=200, .parallel=run.par)

tab_fig4 <- mutate(runs_fig4, ms=bench.time/1000000) %>%
  dcast(N+k+T+bench.rep+ncolors~bench.expr, value.var="ms")  %>%
  mutate(nvars=N*k+k, hess_df=hess/df) %>%
  gather(stat, ms, c(f,df,hess,colors,setup,hess_df)) %>%
  group_by(N, k, T, stat, nvars) %>%
  summarize(mean=mean(ms))


D2 <- filter(data.frame(tab_fig4), stat %in% c("f", "df", "hess",
                                               "colors", "setup","hess_df"))
D2$stat <- plyr::revalue(D2$stat, c("f"="Function", "df"="Gradient", "hess"="Hessian",
                              "colors"="Partitioning",
                              "setup"="Initialization",
                              "hess_df"="Hessian/Gradient"))

D2$stat <- factor(D2$stat, levels=c("Function","Gradient","Hessian",
                                    "Partitioning","Initialization",
                                    "Hessian/Gradient"))

theme_set(theme_bw())
fig4 <- ggplot(D2, aes(x=N,y=mean, color=as.factor(k), linetype=as.factor(k))) %>%
  + geom_line(size=.4) %>%
  + scale_x_continuous("Number of heterogeneous units") %>%
  + scale_y_continuous("Computation time (milliseconds)") %>%
  + guides(color=guide_legend("k"), linetype=guide_legend("k")) %>%
  + facet_wrap(~stat, scales="free") %>%
  + theme(text=element_text(size=8), legend.position="right")


## Replicate Table 4

cases_tab4 <- expand.grid(k=c(2,3,4),
                     N=c(9, 12, 15))
runs_tab4 <- plyr::adply(cases_tab4, 1, run_test_tab4, reps=20, .parallel=run.par)

tab4 <-  mutate(runs_tab4, ms=bench.time/1000000) %>%
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


