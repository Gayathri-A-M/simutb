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
                              nsim = 20000, interval = c(0, 10),
                              ...,
                              verbose = 0,
                              seed = 1000) {

    ## optimization target function
    f_opt <- function(l) {
        t_pfs <- stb_surv_join_simu(nsim, hazard_os, l, rho, seed = seed)
        t_pfs <- t_pfs[, "t_pfs"]

        ## squared loss
        rst <- (median(t_pfs) - median_pfs)^2

        if (verbose > 0)
            cat("median pfs = ",
                median(t_pfs),
                "\n")
        rst
    }

    ## hazard os
    hazard_os <- stb_tl_hazard(median_surv = median_os)

    ## optimization
    rst_optim   <- optimize(f_opt, interval = interval, ...)
    hazard_prog <- rst_optim$minimum

    ## kendall's tau
    t_sim   <- stb_surv_join_simu(nsim, hazard_os, hazard_prog, rho)
    est_cor <- cor.test(t_sim[, "t_pfs"],
                        t_sim[, "t_os"],
                        method = "kendall")$estimate

    rst <- list(median_os   = median_os,
                median_pfs  = median_pfs,
                rho         = rho,
                hazard_os   = hazard_os,
                hazard_prog = hazard_prog,
                kendall     = est_cor)

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
stb_surv_join_simu <- function(n, hazard_os, hazard_prog, rho, seed = NULL) {

    ## inverse cdf for exponential
    f_exp <- function(c, lambda) {
        - log(1 - c) / lambda
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## random bivariate normal samples
    rnd_smp_1 <- rnorm(n)
    rnd_smp_2 <- rho * rnd_smp_1 + sqrt(1 - rho^2) * rnorm(n)

    ## random cdf
    rnd_cdf_1  <- pnorm(rnd_smp_1)
    rnd_cdf_2  <- pnorm(rnd_smp_2)

    ## time to os
    t_os   <- f_exp(rnd_cdf_1, hazard_os)
    t_prog <- f_exp(rnd_cdf_2, hazard_prog)
    t_pfs  <- apply(cbind(t_prog, t_os), 1, min)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    cbind(t_os = t_os, t_pfs = t_pfs, t_prog = t_prog)
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
                                   ...) {

    dta_enroll <- stb_tl_simu_enroll(n,
                                     enroll_dur_mth,
                                     mth_to_days = 30.4,
                                     date_bos    = date_bos,
                                     ...)

    dta_event  <- stb_surv_join_simu(n,
                                     simu_par$hazard_os,
                                     simu_par$hazard_prog,
                                     simu_par$rho) * mth_to_days

    dta_censor <- stb_tl_rexp(n,
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
                                        primary_event = "os",
                                        seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    data_full <- stb_surv_join_trial_simu(...)

    ## interim analysis
    rst <- NULL
    for (i in seq_len(length(info_frac))) {
        data_interim <- stb_tl_interim_data_2arm(data_full,
                                                 total_events,
                                                 info_frac[i],
                                                 event = primary_event)

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
stb_surv_join_summary <- function(data_interim) {

    f_rej <- function(dat, var_rej) {
        n_reps <- max(dat$rep)
        dat %>%
            filter((!!sym(var_rej)) == 1) %>%
            group_by(rep) %>%
            arrange(inx) %>%
            slice(n = 1) %>%
            ungroup() %>%
            group_by(inx, info_frac) %>%
            summarize(N_Rej = n()) %>%
            mutate(Rej = N_Rej / n_reps) %>%
            ungroup() %>%
            mutate(CumuRej = cumsum(Rej))
    }


    rst_events <- data_interim %>%
        group_by(inx, info_frac) %>%
        summarize(nevent_os = mean(nevent_os),
                  nevent_pfs = mean(nevent_pfs))

    rst_os  <- f_rej(data_interim, "rej_os")
    rst_pfs <- f_rej(data_interim, "rej_pfs")

    rst_os$event  <- "os"
    rst_pfs$event <- "pfs"

    list(rejection = rbind(rst_os, rst_pfs),
         nevents   = rst_events)
}
