## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating joint PFS and OS time. All
## functions start with "stb_surv_join"
##
##
##
## -----------------------------------------------------------------------------


#' Calculate hazard and coefficient for simulation
#'
#' Calculate hazard using copula correlation
#'
#' @param nsim number of simulated subjects for numerical calculation
#' @param rho correlation in copula model
#'
#' @export
#'
stb_surv_join_par <- function(median_os, median_pfs, rho,
                              nsim        = 20000,
                              interval_ub = 10,
                              verbose     = 0,
                              method      = c("given.pfs", "given.os"),
                              seed        = 1000,
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

        ## squared loss
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
#' @export
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
#' @export
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


#' Simulate an arm
#'
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
stb_surv_join_arm_simu <- function(n,
                                   simu_par,
                                   enroll_dur_mth,
                                   annual_drop = 0.05,
                                   mth_to_days = 30.4,
                                   date_bos    = as.Date("2022-1-1"),
                                   method      = c("given.pfs", "given.os"),
                                   ...) {

    ## simulation method
    method       <- match.arg(method)
    f_simu_event <- switch(method,
                           given.pfs = stb_surv_join_simu_given_pfs,
                           given.os  = stb_surv_join_simu_given_os)


    dta_enroll <- stb_tl_simu_enroll(n,
                                     enroll_dur_mth,
                                     mth_to_days = 30.4,
                                     date_bos    = date_bos,
                                     ...)

    dta_event  <- f_simu_event(hazard_prog = simu_par$hazard_prog,
                               hazard_os   = simu_par$hazard_os,
                               hazard_pfs  = simu_par$hazard_pfs,
                               rho         = simu_par$rho,
                               nsim        = n, ...) * mth_to_days

    dta_censor <- stb_tl_rexp(n,
                              median_mth  = NULL,
                              hazard      = NULL,
                              annual_drop = annual_drop,
                              mth_to_days = mth_to_days)

    days <- stb_tl_pfs_os(floor(dta_event[, "t_prog"]),
                          floor(dta_event[, "t_os"]),
                          floor(dta_censor))

    cbind(dta_enroll, days) %>%
        data.frame() %>%
        mutate(date_prog   = date_enroll + day_prog,
               date_dth    = date_enroll + day_dth,
               date_censor = date_enroll + day_censor,
               date_pfs    = date_enroll + day_pfs,
               date_os     = date_enroll + day_os,
               date_event  = date_enroll + day_event)
}

#' Simulate a trial
#'
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
stb_surv_join_trial_simu <- function(ntrt,
                                     par_trt,
                                     nctl,
                                     par_ctl,
                                     enroll_dur_mth,
                                     annual_drop = 0.05,
                                     ...) {

    dta_trt <- stb_surv_join_arm_simu(n = ntrt, par_trt, enroll_dur_mth,
                                      annual_drop, ...)

    dta_ctl <- stb_surv_join_arm_simu(n = nctl, par_ctl, enroll_dur_mth,
                                      annual_drop, ...)

    dta_trt$arm <- "1"
    dta_ctl$arm <- "0"

    ## return
    rbind(dta_trt, dta_ctl) %>%
        arrange(date_dth)
}

#' Plot a simulated data
#'
#'
#' @export
#'
stb_surv_join_trial_simu_plot <- function(data) {
    data$tmp <- 1

    fit      <- survfit(Surv(day_prog, tmp) ~ arm, data = data)
    plt_prog <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_dth, tmp) ~ arm, data = data)
    plt_dth  <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_pfs, status_pfs) ~ arm, data = data)
    plt_pfs  <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_os, status_os) ~ arm, data = data)
    plt_os  <- ggsurvplot(fit, data = data)

    list(Death       = plt_dth$plot,
         Progression = plt_prog$plot,
         PFS         = plt_pfs$plot,
         OS          = plt_os$plot)
}


#' Simulate a trial with interim
#'
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
stb_surv_join_trial_interim <- function(...,
                                        total_events,
                                        info_frac,
                                        primary = "os",
                                        seed    = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    data_full <- stb_surv_join_trial_simu(...)

    ## check
    ## data_full %>%
    ##     group_by(arm) %>%
    ##     summarize(m_pfs = median(day_pfs) / 30.4,
    ##               m_os  = median(day_os) / 30.4)

    ## interim analysis
    rst <- NULL
    for (i in seq_len(length(info_frac))) {
        data_interim <- stb_tl_interim_data_2arm(data_full,
                                                 total_events,
                                                 info_frac[i],
                                                 event = primary)

        test_os  <- stb_tl_surv_logrank(data = data_interim,
                                        fml = "Surv(day_os, status_os) ~ arm")

        test_pfs <- stb_tl_surv_logrank(data = data_interim,
                                        fml = "Surv(day_pfs, status_pfs) ~ arm")

        rst <- rbind(rst,
                     c(inx        = i,
                       info_frac  = info_frac[i],
                       pval_os    = unname(test_os["pvalue"]),
                       zscore_os  = unname(test_os["zscore"]),
                       nevent_os  = unname(test_os["nevent"]),
                       hr_os      = unname(test_os["hr"]),
                       pval_pfs   = unname(test_pfs["pvalue"]),
                       zscore_pfs = unname(test_pfs["zscore"]),
                       nevent_pfs = unname(test_pfs["nevent"]),
                       hr_pfs     = unname(test_pfs["hr"]))
                     )
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rownames(rst) <- NULL
    data.frame(rst)
}

#' Simulate a trial with interim
#'
#' @param pval_bounds one-sided p-value boundaries. First row for primary and
#'     second row for secondary
#'
#' @export
#'
stb_surv_join_trial_hiertest <- function(interim_rst, pval_bounds) {
    stopifnot(nrow(pval_bounds) == nrow(interim_rst))

    ps_e <- colnames(pval_bounds)
    if (is.null(ps_e)) {
        ps_e <- c("os", "pfs")
    }

    for (i in seq_len(ncol(pval_bounds))) {
        cur_pval <- interim_rst[,
                                paste("pval_", ps_e[i], sep = "")]
        cur_rej  <- cur_pval <= pval_bounds[, i]

        if (i > 1) {
            cur_rej <- cur_rej & last_rej
        }

        interim_rst[, paste("rej_", ps_e[i], sep = "")] <- cur_rej
        last_rej <- cur_rej
    }

    interim_rst
}


#' Summarize simulation results
#'
#' @export
#'
stb_surv_join_summary <- function(data_interim, primary_secondary) {

    f_rej <- function(dat) {
        n_reps <- max(dat$rep)

        ps  <- paste("rej_", primary_secondary, sep = "")
        dat <- dat %>%
            filter((!!sym(ps[1])) == 1) %>%
            group_by(rep) %>%
            arrange(inx) %>%
            slice(n = 1)

        rej_pri <- dat %>%
            group_by(inx, info_frac) %>%
            summarize(N_Rej = n()) %>%
            mutate(Rej = N_Rej / n_reps) %>%
            ungroup() %>%
            mutate(CumuRej = cumsum(Rej))

        rej_sec <- dat %>%
            filter((!!sym(ps[2])) == 1) %>%
            group_by(inx, info_frac) %>%
            summarize(N_Rej = n()) %>%
            mutate(Rej = N_Rej / n_reps) %>%
            ungroup() %>%
            mutate(CumuRej = cumsum(Rej))

        rej_pri$event <- primary_secondary[1]
        rej_sec$event <- primary_secondary[2]

        rbind(rej_pri, rej_sec)
    }


    ## events
    rst_events <- data_interim %>%
        group_by(inx, info_frac) %>%
        summarize(nevent_os  = mean(nevent_os),
                  nevent_pfs = mean(nevent_pfs))

    ## marginal rejection
    rst_rej_marginal <- data_interim %>%
        gather(type, rej, rej_os, rej_pfs) %>%
        group_by(inx, info_frac, type) %>%
        summarize(rej = mean(rej))

    ## rejection probability
    rst_rejection <- f_rej(data_interim)

    ## hazard ratio
    mean_hr <- data_interim %>%
        gather(type, hr, hr_os, hr_pfs) %>%
        group_by(type, inx, info_frac) %>%
        summarize(hr = mean(hr))

    ## effect size
    mean_zscore <- data_interim %>%
        gather(type, zscore, zscore_os, zscore_pfs) %>%
        group_by(type, inx, info_frac) %>%
        summarize(zscore = mean(zscore))

    ## correlation matrix
    n_ana <- max(data_interim$inx)
    z_os <- matrix(data_interim$zscore_os,
                   ncol  = n_ana,
                   byrow =  TRUE)
    z_pfs <- matrix(data_interim$zscore_pfs,
                    ncol  = n_ana,
                    byrow =  TRUE)

    cor_matrix <- cor(cbind(z_os, z_pfs))
    rownames(cor_matrix) <- colnames(cor_matrix) <-
        c(paste("os_",  seq_len(n_ana), sep = ""),
          paste("pfs_", seq_len(n_ana), sep = ""))

    ## correlation of z_primary and z_secondary
    cor_z <- cor_matrix[(n_ana + 1) : (2 * n_ana), 1 : n_ana]
    cor_z <- mean(diag(cor_z))

    ## return
    list(rej_marginal  = rst_rej_marginal,
         rejection     = rst_rejection,
         nevents       = rst_events,
         hr_mean       = mean_hr,
         zscore_mean   = mean_zscore,
         zscore_cormat = cor_matrix,
         zscore_cor    = cor_z)
}
