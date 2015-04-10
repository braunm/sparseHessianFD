
color.list2vec <- function(W, n) {

    k <- length(W)
    res <- rep(0,n)
    for (i in 1:k) {
        res[as.integer(W[[i]])] <- as.integer(i)
    }
    return(res)
}


color.cols <- function(rows, cols) {

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
            Sr <- as.vector(S[r,,drop=FALSE] %*% S)
            nei <- trackS & Sr
            deg[nei] <- 0
            trackS[nei] <- FALSE
            W[[k]] <- c(W[[k]], r)
        }
        G[,W[[k]]] <- FALSE
        trackG[W[[k]]] <- FALSE
    }
    return(W)
}



color.cols.C <- function(rows, cols, nvars) {
    diagElement <- rows == cols
    rows2 <- c(rows, cols[!diagElement]) - 1
    cols2 <- c(cols, rows[!diagElement]) - 1
 ##   rows2 <- rows-1
 ##   cols2 <- cols-1
    W <- color(as.integer(rows2), as.integer(cols2), nvars)
    return(W)
}


