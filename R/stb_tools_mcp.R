
#' Hochberg multiplicity control
#'
#' @export
stb_multi_hochberg <- function(pvals_0, alpha = 0.05, ...) {
    p_inx <- order(pvals_0, decreasing = TRUE)
    pvals <- pvals_0[p_inx]

    rej   <- rep(1, length(pvals))
    for (i in seq_len(length(pvals))) {
        if (pvals[i] < alpha / i)
            break

        rej[i] <- 0
    }

    rej[order(p_inx)]
}

#' Holmes multiplicity control
#'
#' @export
stb_multi_holms <- function(pvals, alpha = 0.05, ...) {
    p_inx <- order(pvals)
    pvals <- pvals[p_inx]
    k     <- length(pvals)

    rej   <- rep(0, k)
    for (i in seq_len(k)) {
        if (pvals[i] > alpha / (k - i + 1))
            break

        rej[i] <- 1
    }

    rej[order(p_inx)]
}

#' Hierarchical
#'
#' @export
stb_multi_hierarchi <- function(pvals, alpha  = 0.05, p_inx = NULL,  ...) {
    k <- length(pvals)
    if (is.null(p_inx)) {
        p_inx <- seq_len(k)
    }

    pvals <- pvals[p_inx]
    rej   <- rep(0, k)
    for (i in seq_len(k)) {
        if (pvals[i] > alpha)
            break

        rej[i] <- 1
    }

    rej[order(p_inx)]
}
