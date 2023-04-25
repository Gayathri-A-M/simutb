## -----------------------------------------------------------------------------
##
## This file contains R functions for distributions
##
##
##
## -----------------------------------------------------------------------------


#' Optimize to get parameters
#'
#'
internal_opt_para <- function(target_qs, para, f_qs, ...) {

    fn <- function(para) {
        cur_val <- do.call(f_qs, c(list(p = ps), para))
        sum((cur_val - val)^2)
    }

    ps  <- as.numeric(names(target_qs))
    val <- target_qs
    rst <- optim(para, fn, ...)
    rst
}


#' Get beta parameters by quantiles
#'
#'  @export
#'
stb_dist_beta_para_by_qs <- function(target_qs = c("0.1" = 0.5,
                                                   "0.5" = 0.7),
                                     para      = c(shape1 = 1, shape2 = 1),
                                     lower     = c(0, 0),
                                     method    = "L-BFGS-B",
                                     ...) {

    rst <- internal_opt_para(target_qs,
                             para   = para,
                             f_qs   = qbeta,
                             lower  = lower,
                             method = method,
                             ...)

    if (0 != rst$convergence) {
        cat("Convergence was not achieved. \n")
        rst <- list()
    } else {
        rst <- list(shape1 = unname(rst$par[1]),
                    shape2 = unname(rst$par[2]))
    }

    rst
}

#' Get normal parameters by quantiles
#'
#'  @export
#'
stb_dist_normal_para_by_qs <- function(target_qs = c("0.1" = 0.5,
                                                     "0.5" = 0.7),
                                       para      = c(mean = 1, sd = 1),
                                       lower     = c(-Inf, 0),
                                       method    = "L-BFGS-B",
                                       ...) {

    rst <- internal_opt_para(target_qs,
                             para   = para,
                             f_qs   = qnorm,
                             lower  = lower,
                             method = method,
                             ...)

    if (0 != rst$convergence) {
        cat("Convergence was not achieved. \n")
        rst <- list()
    } else {
        rst <- list(mean = unname(rst$par[1]),
                    sd   = unname(rst$par[2]))
    }

    rst
}


#' Get log-normal parameters by quantiles
#'
#'  @export
#'
stb_dist_lognormal_para_by_qs <- function(target_qs = c("0.1" = 0.5,
                                                        "0.5" = 0.7),
                                          ...) {
    target_qs <- log(target_qs)
    par_norm  <- stb_dist_normal_para_by_qs(target_qs, ...)
    par_logn  <- tl_lognorm_par(par_norm$mean,
                                par_norm$sd)

    list(mean       = par_logn$mean,
         sd         = par_logn$sd,
         norm_mean  = par_norm$mean,
         norm_sd    = par_norm$sd)
}
