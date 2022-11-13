#' Calculate hazard and coefficient for simulation
#'
#' Calculate hazard using copula correlation
#'
#' @param nsim number of simulated subjects for numerical calculation
#' @param rho correlation in copula model
#'
#'
stb_surv_join_par_v0 <- function(median_os, median_pfs, rho,
                              nsim        = 20000,
                              interval_ub = 10,
                              verbose     = 0,
                              method      = c("given.pfs", "given.os"),
                              seed        = NULL,
                              nlarge      = 50000,
                              ...) {

    ## optimization target function
    f_opt_given_os <- function(l) {
        t_sim <- stb_surv_join_simu_given_os(l,
                                             hazard_os, hazard_pfs, rho,
                                             rnd_cdf = rnd_cdf)
        t_pfs <- t_sim[, "t_pfs"]

        ## squared loss
        rst <- (median(t_pfs) - median_pfs)^2

        if (verbose > 0)
            cat("median pfs = ",
                median(t_pfs),
                "/", median_pfs,
                "\n")
        rst
    }

    ## optimization target function
    f_opt_given_pfs <- function(l) {
        t_sim <- stb_surv_join_simu_given_pfs(l,
                                              hazard_os, hazard_pfs, rho,
                                              rnd_cdf = rnd_cdf,
                                              t_pfs = t_pfs, ...)

        t_os <- t_sim[, "t_os"]

        ## loss
        rst <- median(t_os) - median_os

        if (verbose > 0)
            cat("median os = ",
                median(t_os),
                "/", median_os,
                ", hazard_prog = ",
                l,
                "\n")
        rst
    }

    ## optimization method
    method <- match.arg(method)

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## hazard
    hazard_os  <- stb_tl_hazard(median_surv = median_os)
    hazard_pfs <- stb_tl_hazard(median_surv = median_pfs)
    t_pfs      <- rexp(nsim, hazard_pfs)

    ## optimization
    rnd_cdf     <- stb_tl_binorm(nlarge, rho)$rnd_cdf
    f_opt       <- switch(method,
                          given.pfs = f_opt_given_pfs,
                          given.os  = f_opt_given_os)

    rst_optim   <- uniroot(f_opt,
                           interval = c(0, interval_ub))
    hazard_prog <- rst_optim$root

    ## kendall's tau
    if ("given.os" == method) {
        t_sim   <- stb_surv_join_simu_given_os(hazard_prog,
                                               hazard_os,
                                               hazard_pfs,
                                               rho,
                                               nsim = nsim)
    } else {
        t_sim   <- stb_surv_join_simu_given_pfs(hazard_prog,
                                                hazard_os,
                                                hazard_pfs,
                                                rho,
                                                nsim = nsim,
                                                ...)
    }
    est_cor <- cor.test(t_sim[, "t_pfs"],
                        t_sim[, "t_os"],
                        method = "kendall")$estimate

    if (!is.null(seed))
        set.seed(old_seed)

    rst <- list(median_os   = median_os,
                median_pfs  = median_pfs,
                rho         = rho,
                hazard_os   = hazard_os,
                hazard_pfs  = hazard_pfs,
                hazard_prog = hazard_prog,
                kendall     = est_cor,
                optim_raw   = rst_optim)

    class(rst) <- "join_surv_par"

    rst
}


#' Simulate OS and PFS by Copula and Exponential
#'
#' @param n sample size
#' @param hazard_os overall survival hazard
#' @param hazard_prog  progression hazard
#' @param rho correlation between OS and progression in gaussian copula
#' @param rnd_cdf exisiting samples
#'
#'
stb_surv_join_simu_given_os <- function(hazard_prog, hazard_os, hazard_pfs, rho,
                                        rnd_cdf = NULL, nsim,
                                        seed = NULL, ...) {

    ## inverse cdf for exponential
    f_exp <- function(c, lambda) {
        - log(1 - c) / lambda
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    if (is.null(rnd_cdf)) {
        rnd_cdf <- stb_tl_binorm(nsim, rho)$rnd_cdf
    }

    ## time to os. rnd_cdf may have more rows than n
    t_os   <- f_exp(rnd_cdf[1:nsim, 1], hazard_os)
    t_prog <- f_exp(rnd_cdf[1:nsim, 2], hazard_prog)
    t_pfs  <- apply(cbind(t_prog, t_os), 1, min)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    cbind(t_os = t_os, t_pfs = t_pfs, t_prog = t_prog)
}

#' Simulate OS and PFS by Copula and Exponential
#'
#' @param n sample size
#' @param hazard_os overall survival hazard
#' @param hazard_prog  progression hazard
#' @param rho correlation between OS and progression in gaussian copula
#' @param rnd_cdf exisiting samples
#'
#'
stb_surv_join_simu_given_pfs <- function(hazard_prog, hazard_os, hazard_pfs,
                                         rho,
                                         rnd_cdf = NULL, t_pfs = NULL,
                                         nsim = NULL,
                                         h = 3, seed = NULL, nlarge = 20000,
                                         n_core = 1) {

    f_os_given_pfs <- function(t, smp_os, smp_prog,
                               hazard_prog, hazard_os, hazard_pfs,
                               h) {

        ## prob t- h < pfs < t + h
        p_denom <- exp(-hazard_pfs * (t - h)) - exp(-hazard_pfs * (t + h))

        ## prob t- h < os  < t + h & prog > os
        inx_num <- which(
            smp_os   > t - h &
            smp_os   < t + h &
            smp_prog > smp_os)

        p_num <- length(inx_num) / length(smp_os)

        ## is current pfs time t os?
        p_os_lt_prog <- min(p_num / p_denom, 1)
        ind          <- rbinom(1, 1, p_os_lt_prog)

        ## os time given pfs = t
        if (1 == ind) {
            ## os < prog
            rst_os   <- t
            rst_prog <- t + rexp(1, hazard_prog)
        } else {
            ## os > prog
            rst_os   <- t + rexp(1, hazard_os)
            rst_prog <- t
        }

        c(rst_os, rst_prog)
    }


    ## inverse cdf for exponential
    f_exp <- function(c, lambda) {
        - log(1 - c) / lambda
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    if (is.null(rnd_cdf)) {
        rnd_cdf <- stb_tl_binorm(nlarge, rho)$rnd_cdf
    }

    if (is.null(t_pfs)) {
        t_pfs <- rexp(nsim, hazard_pfs)
    }

    ## time to os
    large_smp_os   <- f_exp(rnd_cdf[, 1], hazard_os)
    large_smp_prog <- f_exp(rnd_cdf[, 2], hazard_prog)

    ## conditional os
    t_cond <- parallel::mclapply(t_pfs,
                                 function(x) {
                                     f_os_given_pfs(
                                         x,
                                         smp_os      = large_smp_os,
                                         smp_prog    = large_smp_prog,
                                         hazard_prog = hazard_prog,
                                         hazard_os   = hazard_os,
                                         hazard_pfs  = hazard_pfs,
                                         h           = h)
                                 },
                                 mc.cores = n_core)

    t_cond <- simplify2array(t_cond)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    cbind(t_os   = t_cond[1, ],
          t_pfs  = t_pfs,
          t_prog = t_cond[2, ])
}
