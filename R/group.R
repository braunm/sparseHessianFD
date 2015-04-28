
#' @param M full pattern matrix
get_groups <- function(M) {

    L <- tril(as(M,"nMatrix"))
    G <- Matrix::crossprod(L)

    nvars <- NROW(M)
    forb <- vector("list", length=nvars) ## forbidden colors
    colors <- rep(0, nvars)

    for (i in 1:nvars) {
        ## color = smallest integer not forbidden for i
        if (is.null(forb[[i]])) {
            colors[i] <- 1
        } else {
            colors[i] <- min(seq(1:(max(forb[[i]])+1))[-forb[[i]]])
        }

        for (j in which(G[i,])) forb[[j]] <- c(forb[[j]], colors[i])
    }

    return(colors)
}



