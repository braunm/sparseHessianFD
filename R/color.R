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

color.list2vec <- function(W, n) {

    k <- length(W)
    res <- rep(0,n)
    for (i in 1:k) {
        res[as.integer(W[[i]])] <- i
    }
    return(res)
}




