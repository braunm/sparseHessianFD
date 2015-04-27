
#' @param M full pattern matrix
get_groups <- function(M) {

    X <- as(M, "nMatrix")
    P <- order(Matrix::rowSums(X), decreasing=TRUE)
    L <- tril(X[P,P])
    tmp <- L %&% L ## intersection graph from Boolean matrix product


    G <-  as(symmpart(tmp), "nMatrix")


    nvars <- NROW(M)
    forb <- vector("list", length=nvars) ## forbidden colors
    colors <- rep(0, nvars)

    for (i in 1:nvars) {
        ## color = smallest integer not forbidden for i
        if (is.null(forb[[i]])) {
            colors[i] <- 0
        } else {
            colors[i] <- min(seq(0:(max(forb[[i]])+1))[-forb[[i]]])
        }
        for (j in which(G[i,])) forb[[j]] <- c(forb[[j]], colors[i])
    }

    return(colors)
}


