rm(list=ls())

library(sparseHessianFD)
library(dplyr)
library(igraph)

color.cols <- function(rows, cols) {
    G <- graph(rbind(rows, cols))
    nV <- vcount(G)
    ## Since igraph renumbers vertices after some are removed,
    ## we need to track the vertices using names.
    G <- set.vertex.attribute(G, "name", value=as.character(1:nV))

    k <- 0
    W <- NULL ##vector("list", length=1)
    while (vcount(G) > 0) {
        k <- k + 1
  ##      cat("\nk = ",k,"\n")
  ##      print(V(G))
        S <- G
        nm <- get.vertex.attribute(S, "name")
        W <- c(W, list(NULL))
        while (vcount(S) > 0) {
            deg <- degree(S, loops=FALSE, mode="out")
            r <- nm[which.max(deg)]
            ##       cat("\tmax.degree.vertex = ",r,"\n")
            W[[k]] <- c(W[[k]], r)
            nei <- nm[neighborhood(S, 2, r)[[1]]]
            ##       cat("\tneighborhood = ",nei,"\n")
            S <- delete.vertices(S, nei)
            nm <- get.vertex.attribute(S, "name")
        }
        G <- G - W[[k]]
  ##      cat("\tW[[",k,"]] = ", W[[k]], "\n")
    }
    return(W)
}

## Test matrices H and L
b <- 3
n <- 22
p <- 2

H <- Matrix(F, nrow=n*b+p, ncol=n*b+p)
for (j in 1:b) {
    for (i in 1:n) {
        H[seq(i,n*b,by=n),(j-1)*n+i] <- T
        H[(n*b+1):(n*b+p), (j-1)*n+i] <- T
    }
}
H[,(n*b+1):(n*b+p)] <- T


L <- bdiag(replicate(n, matrix(T,b,b), simplify=FALSE)) %>%
  rBind(matrix(T,p,n*b)) %>%
  cBind(matrix(T, b*n+p,p))


testmat <- L
hs <- Matrix.to.Coord(tril(testmat))

W <- color.cols(hs$rows, hs$cols)



