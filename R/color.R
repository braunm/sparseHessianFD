## color.cols2 <- function(rows, cols) {
##     G <- graph(rbind(rows, cols))
##     nV <- vcount(G)
##     ## Since igraph renumbers vertices after some are removed,
##     ## we need to track the vertices using names.
##     G <- set.vertex.attribute(G, "name", value=as.character(1:nV))

##     k <- 0
##     W <- NULL ##vector("list", length=1)
##     while (vcount(G) > 0) {
##         k <- k + 1
##   ##      cat("\nk = ",k,"\n")
##   ##      print(V(G))
##         S <- G
##         nm <- get.vertex.attribute(S, "name")
##         W <- c(W, list(NULL))
##         while (vcount(S) > 0) {
##             deg <- degree(S, loops=FALSE, mode="out")
##             r <- nm[which.max(deg)]
##             ##       cat("\tmax.degree.vertex = ",r,"\n")
##             W[[k]] <- c(W[[k]], r)
##             nm.idx <- neighborhood(S, 2, r)
##             nei <- nm[nm.idx[[1]]]
##             ##       cat("\tneighborhood = ",nei,"\n")
##             S <- delete.vertices(S, nei)
##             nm <- get.vertex.attribute(S, "name")
##         }
##         G <- G - W[[k]]
##   ##      cat("\tW[[",k,"]] = ", W[[k]], "\n")
##     }
##     for (i in 1:length(W)) {
##         W[[i]] <- as.integer(W[[i]])
##     }
##     return(W)
## }


color.list2vec <- function(W, n) {

    k <- length(W)
    res <- rep(0,n)
    for (i in 1:k) {
        res[as.integer(W[[i]])] <- as.integer(i)
    }
    return(res)
}


color.cols1 <- function(rows, cols) {

    A <- sparseMatrix(i=rows, j=cols) ## lower triangle
    stopifnot(Matrix::isTriangular(A, upper=FALSE))
    G <- as(tril(A) + Matrix::t(tril(A, k=-1)), "nMatrix") ## symmetric
    stopifnot(Matrix::isSymmetric(G))

    nV <- NROW(G)
    dimnames(G) <- list(1:nV, 1:nV)

    k <- 0
    W <- NULL
    while (NROW(G) > 0) {
        k <- k + 1
        S <- G
        nm <- rownames(S)
        W <- c(W, list(NULL))
        while (NROW(S) > 1) {
            deg <- Matrix::rowSums(S)-1
            r <- nm[which.max(deg)]
            W[[k]] <- c(W[[k]], r)
            S2 <- as.logical(S[r,,drop=FALSE] %*% S)
            nei <- nm[S2]
            d <- which(nm %in% nei)
            S <- S[-d, -d, drop=FALSE]
            nm <- nm[-d]
        }
        if (NROW(S) <= 1) {
            W[[k]] <- c(W[[k]], nm[1])
            S <- NULL
            nm <- NULL
        }
        d <- which(rownames(G) %in% W[[k]])
        G <- G[-d, -d, drop=FALSE]
    }
    for (i in 1:length(W)) {
        W[[i]] <- sort(as.integer(W[[i]]), decreasing=FALSE)
    }
    return(W)
}

color.cols2 <- function(rows, cols) {

    A <- sparseMatrix(i=rows, j=cols) ## lower triangle
    stopifnot(Matrix::isTriangular(A, upper=FALSE))
    G <- as(tril(A) + Matrix::t(tril(A, k=-1)), "nMatrix") ## symmetric
    stopifnot(Matrix::isSymmetric(G))

    nV <- NROW(G)
    k <- 0
    W <- NULL
    trackG <- rep(TRUE, nV)
    while (any(trackG)) {
        k <- k + 1
        deg <- Matrix::rowSums(G) * trackG
        W <- c(W, list(NULL))
        S <- G
        trackS <- trackG
        while (any(trackS)) {
            r <- which.max(deg)
            W[[k]] <- c(W[[k]], r)
            Sr <- S[r,,drop=FALSE] %*% S
            nei <- trackS & as.logical(Sr)
            deg[nei] <- 0
            trackS[nei] <- FALSE
        }
        G[W[[k]],] <- FALSE
        G[,W[[k]]] <- FALSE
        trackG[W[[k]]] <- FALSE
    }
    for (i in 1:length(W)) {
        W[[i]] <- sort(as.integer(W[[i]]), decreasing=FALSE)
    }
    return(W)
}


color.cols <- color.cols2

