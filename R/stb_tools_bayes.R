## -----------------------------------------------------------------------------
##
## This file contains R tool functions for survival analysis.
##
## All functions start with "stb_tl_".
##
##
##
## -----------------------------------------------------------------------------




#' Get posterior distribution
#'
#'
#' @export
#'
stb_tl_bayes_beta_post <- function(obs_y  = 1,
                                   obs_n  = 10,
                                   pri_ab = c(1, 1),
                                   x      = seq(0, 1, by = 0.01),
                                   n_post = NULL) {
    pos_a   <- pri_ab[1] + obs_y
    pos_b   <- pri_ab[2] + obs_n - obs_y

    rst_den <- NULL
    if (!is.null(x)) {
        rst <- sapply(x, function(v) {
            c(dbeta(v, pos_a, pos_b),
              1 - pbeta(v, pos_a, pos_b))
        })

        rst_den <- cbind(x     = x,
                         y_pdf = rst[1, ],
                         y_cdf = rst[2, ])
    }

    post_smp <- NULL
    if (!is.null(n_post)) {
        post_smp <- rbeta(n_post, pos_a, pos_b)
    }


    list(pri_ab   = pri_ab,
         pos_ab   = c(pos_a, pos_b),
         density  = rst_den,
         post_smp = post_smp)
}


#' Elicit prior
#'
#'
#' @export
#'
stb_tl_bayes_elicit <- function(x_range = c(0.8, 0.9),
                                x_prob  = 0.9,
                                ub      = 1000) {
    fp <- function(b) {
        a     <- b * mu / (1 - mu)
        cur_p <- pbeta(x_range[2], a, b) - pbeta(x_range[1], a, b)
        cur_p - x_prob
    }

    mu <- mean(x_range)
    b  <- uniroot(fp, c(0, ub))$root

    c(pri_a = b * mu / (1 - mu),
      pri_b = b)
}


#' Generate Binomial Outcome
#'
#'
#' @export
#'
stb_tl_bayes_gen <- function(n, p, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)


    y <- rbinom(1, n, p)

    if (!is.null(seed))
        set.seed(old_seed)

    cbind(n = n, p = p, y = y)
}


#' Bayesian decision making
#'
#' @export
#'
stb_tl_bayes_decision <- function(post_smp,
                                  decision_h0,
                                  decision_thresh,
                                  decision_gl = c("greater.than",
                                                  "less.than")) {

    decision_gl   <- match.arg(decision_gl)
    decision_prob <- switch(decision_gl,
                            greater.than = mean(post_smp > decision_h0),
                            less.than    = mean(post_smp < decision_h0))

    c(success     = decision_prob > decision_thresh,
      probability = decision_prob)
}
