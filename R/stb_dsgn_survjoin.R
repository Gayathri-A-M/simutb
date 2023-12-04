## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating joint PFS and OS time. All
## functions start with "stb_surv_join"
##
##
##
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
survjoin_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Oncology trial joint PFS and OS as the primary and \n")
    cat("    secondary endpoints. \n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 500)\n")
    cat("    ctl_ratio:      randomization proportion for control \n")
    cat("                    (default 0.5)\n")
    cat("    info_frac:      information fractions for interim analysis \n")
    cat("                    (default c(0.3, 0.7, 1))\n")
    cat("    primary:        primary endpoint, either OS (os) or PFS (pfs) \n")
    cat("                    (default os)\n")
    cat("    target_primary: target number of primary events for the \n")
    cat("                    final analysis (default 100)\n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    btype_prmary:   group sequential boundary type for the \n")
    cat("                    primary endpoint  (default asOF)\n")
    cat("    bound_second:   group sequential boundary for the secondary\n")
    cat("                    endpoint  (default same as primary endpoint)\n")
    cat("    ctl_median_os:  control median OS in months (default 15)\n")
    cat("    ctl_median_pfs: control median PFS in months (default 15)\n")
    cat("    hr_os:          hazard ratio (treatment vs. control) for OS\n")
    cat("                    (default 0.7)\n")
    cat("    hr_pfs:         hazard ratio (treatment vs. control) for PFS\n")
    cat("                    (default 1)\n")
    cat("    rho_ctl:        control group correlation between PFS and OS\n")
    cat("                    test statistic (default 0.8) \n")
    cat("    rho_trt:        treatment group correlation between PFS and OS\n")
    cat("                    test statistic (default the same as control) \n")
    cat("    bound_second:   group sequential boundary for the secondary\n")
    cat("                    endpoint  (default same as primary endpoint)\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    enroll_dur_mth: enrollment months (default 18) \n")
    cat("    date_bos:       begin of study date\n")
    cat("                    (default 2022-1-1) \n")
    cat("    seed:           random seed for calculating the design para\n")
}


#' Set design paramter
#'
#'
survjoin_default_para <- function(x) {
    lst_default <- list(info_frac       = c(0.3, 0.7, 1),
                        primary         = "os",
                        target_primary  = 100,
                        sample_size     = 500,
                        ctl_ratio       = 0.5,
                        annual_drop     = 0.000001,
                        enroll_dur_mth  = 18,
                        alpha           = 0.025,
                        btype_primary   = "asOF",
                        bound_second    = NULL,
                        ctl_median_os   = 15,
                        ctl_median_pfs  = 9,
                        hr_os           = 0.7,
                        hr_pfs          = 1,
                        rho_ctl         = 0.1,
                        rho_trt         = NULL,
                        seed            = NULL)

    do.call(internal_survjoin_dpara, lst_default)
}

#' Get Study Design for Joint Primary and Secondary Endpoints
#'
internal_survjoin_dpara <- function(info_frac       = c(0.3, 0.7, 1),
                                    primary         = c("os", "pfs"),
                                    target_primary  = 100,
                                    sample_size     = 500,
                                    ctl_ratio       = 0.5,
                                    annual_drop     = 0.000001,
                                    enroll_dur_mth  = 18,
                                    alpha           = 0.025,
                                    btype_primary   = "asOF",
                                    bound_second    = NULL,
                                    ctl_median_os   = 15,
                                    ctl_median_pfs  = 9,
                                    hr_os           = 0.7,
                                    hr_pfs          = 1,
                                    rho_ctl         = 0.8,
                                    rho_trt         = NULL,
                                    method          = c("given.pfs",
                                                        "given.os"),
                                    ...,
                                    seed            = NULL) {

    primary <- match.arg(primary)
    method  <- match.arg(method)

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)


    ## hazard and correlation
    par_ctl <- stb_surv_join_par(median_os  = ctl_median_os,
                                 median_pfs = ctl_median_pfs,
                                 rho        = rho_ctl)

    if (is.null(rho_trt))
        rho_trt <- rho_ctl

    par_trt <- stb_surv_join_par(median_os  = ctl_median_os  / hr_os,
                                 median_pfs = ctl_median_pfs / hr_pfs,
                                 rho        = rho_trt)

    ## boundary
    bound_primary <-
        getDesignGroupSequential(sided            = 1,
                                 alpha            = alpha,
                                 informationRates = info_frac,
                                 typeOfDesign     = btype_primary)$stageLevels

    if (is.null(bound_second))
        bound_second  <- rep(alpha, length(bound_primary))

    pval_bounds           <- cbind(bound_primary, bound_second)
    colnames(pval_bounds) <- unique(c(primary, c("os", "pfs")))

    ## sample size
    nctl <- floor(sample_size * ctl_ratio)
    ntrt <- sample_size - nctl

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    as.list(environment())
}


#' Simulate a trial
#'
#'
survjoin_gen_data <- function(ntrt,
                              par_trt,
                              nctl,
                              par_ctl,
                              ...) {

    dta_trt <- stb_surv_join_arm_simu(n = ntrt, par_trt, ...)
    dta_ctl <- stb_surv_join_arm_simu(n = nctl, par_ctl, ...)

    dta_trt$arm <- "1"
    dta_ctl$arm <- "0"

    ## return
    rbind(dta_trt, dta_ctl) %>%
        arrange(date_dth)
}

#' Analyze a trial
#'
#'
#'
survjoin_ana_data <- function(data_full,
                              total_events,
                              info_frac,
                              primary = "os") {

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
                       pval_os    = unname(test_os["pvalue_oneside"]),
                       zscore_os  = unname(test_os["zscore"]),
                       nevent_os  = unname(test_os["nevent"]),
                       hr_os      = unname(test_os["hr"]),
                       pval_pfs   = unname(test_pfs["pvalue_oneside"]),
                       zscore_pfs = unname(test_pfs["zscore"]),
                       nevent_pfs = unname(test_pfs["nevent"]),
                       hr_pfs     = unname(test_pfs["hr"]))
                     )
    }

    ## return
    rownames(rst) <- NULL
    data.frame(rst)
}


#' Plot a simulated data
#'
#'
#'
survjoin_plot_data <- function(data) {
    data$tmp <- 1
    fit      <- survfit(Surv(day_prog, tmp) ~ arm, data = data)
    plt_prog <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_dth, tmp) ~ arm, data = data)
    plt_dth  <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_pfs, status_pfs) ~ arm, data = data)
    plt_pfs  <- ggsurvplot(fit, data = data)

    fit      <- survfit(Surv(day_os, status_os) ~ arm, data = data)
    plt_os   <- ggsurvplot(fit, data = data)

    list(Death       = plt_dth$plot,
         Progression = plt_prog$plot,
         PFS         = plt_pfs$plot,
         OS          = plt_os$plot)
}



#' Simulate a trial with interim
#'
#'
#'
survjoin_hiertest <- function(interim_rst, pval_bounds) {
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
survjoin_simu_summary <- function(data_interim, primary_secondary) {

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
    cor_z <- NULL
    for (i in 1:n_ana) {
        cur_cor <- cor_matrix[paste("os_",  i, sep = ""),
                              paste("pfs_", i, sep = "")]
        cor_z   <- c(cor_z, cur_cor)
    }
    cor_z <- mean(cor_z)

    ## return
    list(rej_marginal  = rst_rej_marginal,
         rejection     = rst_rejection,
         nevents       = rst_events,
         hr_mean       = mean_hr,
         zscore_mean   = mean_zscore,
         zscore_cormat = cor_matrix,
         zscore_cor    = cor_z)
}

#' Summarize key simulation results
#'
#' @export
#'
survjoin_simu_key <- function(rst_summary, lst_design) {
    f_rej <- function(endp) {
        tmp <- rst_summary$rejection %>%
            filter(event == endp)

        if (0 == nrow(tmp)) {
            rst <- 0
        } else {
            rst <- max(tmp$CumuRej)
        }

        rst
    }

    rst_key <- data.frame(
        pri_sec        = paste(colnames(lst_design$pval_bounds),
                               collapse = ","),
        info_frac      = paste(lst_design$info_frac,
                               collapse = ","),
        hr_os          = lst_design$hr_os,
        hr_pfs         = lst_design$hr_pfs,
        sample_size    = lst_design$sample_size,
        target_primary = lst_design$target_primary,
        rho_ctl        = lst_design$rho_ctl,
        rho_trt        = lst_design$rho_trt,
        hz_prog_ctl    = lst_design$par_ctl$hazard_prog,
        hz_prog_trt    = lst_design$par_trt$hazard_prog,
        zscore_cor     = rst_summary$zscore_cor,
        rej_os         = f_rej("os"),
        rej_pfs        = f_rej("pfs"))


    rst_key
}

#' Calculate hazard and coefficient for simulation
#'
#'
#' @export
#'
stb_surv_join_par <- function(median_os, median_pfs, rho) {

    f_lv <- function(lu, hazard_pfs, rho) {
        sqrt(rho^2 / (1 - rho^2) * (hazard_pfs^2 - lu^2))
    }

    f_opt <- function(lu) {
        lp <- hazard_pfs - lu
        lv <- f_lv(lu, hazard_pfs, rho)

        t1 <- lp / (hazard_pfs - lv) * exp(-lv * median_os)
        t2 <- (lv - lu) / (hazard_pfs - lv) * exp(-hazard_pfs * median_os)

        1 - t1 + t2 - 1/2
    }

    stopifnot(median_os >= median_pfs)

    hazard_pfs  <- stb_tl_hazard(median_surv = median_pfs)
    rst_optim   <- uniroot(f_opt, interval = c(0, hazard_pfs))
    hazard_u    <- rst_optim$root
    hazard_prog <- hazard_pfs - hazard_u
    hazard_os   <- f_lv(hazard_u, hazard_pfs, rho)

    rst <- list(median_os   = median_os,
                median_pfs  = median_pfs,
                rho         = rho,
                hazard_os   = hazard_os,
                hazard_pfs  = hazard_pfs,
                hazard_prog = hazard_prog,
                hazard_u    = hazard_u,
                optim_raw   = rst_optim)

    class(rst) <- "PAR_PFS_OS"

    rst
}

#' Simulate OS and PFS by Copula and Exponential
#'
#' @param n sample size
#' @param hazard_os overall survival hazard
#' @param hazard_prog  progression hazard
#' @param hazard_pfs   progression free survival hazard
#'
#' @export
#'
stb_surv_join_simu_events <- function(n, par_pfs_os, seed = NULL) {

    stopifnot("PAR_PFS_OS" %in% class(par_pfs_os))

    if (!is.null(seed))
        old_seed <- set.seed(seed)


    t_prog  <- rexp(n, par_pfs_os$hazard_prog)
    t_u     <- rexp(n, par_pfs_os$hazard_u)
    t_pfs   <- apply(cbind(t_prog, t_u), 1, min)

    ## OS
    t_os     <- t_pfs
    inx      <- which(t_prog == t_pfs)
    if (length(inx) > 0)
        t_os[inx] <- t_os[inx] + rexp(length(inx), par_pfs_os$hazard_os)


    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    data.frame(t_os   = t_os,
               t_pfs  = t_pfs,
               t_prog = t_prog,
               t_u    = t_u)
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
                                   annual_drop = 0,
                                   mth_to_days = 30.4,
                                   date_bos    = as.Date("2022-1-1"),
                                   ...) {

    dta_enroll <- stb_tl_simu_enroll_arc(n,
                                         enroll_dur_mth,
                                         date_bos    = date_bos,
                                         ...)

    dta_censor <- stb_tl_rexp(n,
                              median_mth  = NULL,
                              hazard      = NULL,
                              annual_drop = annual_drop,
                              ...)

    dta_event  <- stb_surv_join_simu_events(n, simu_par)
    days       <- stb_tl_pfs_os(dta_event[, "t_prog"],
                                dta_event[, "t_os"],
                                dta_censor)

    cbind(dta_enroll, days) %>%
        data.frame() %>%
        mutate(date_prog   = date_enroll + day_prog,
               date_dth    = date_enroll + day_dth,
               date_censor = date_enroll + day_censor,
               date_pfs    = date_enroll + day_pfs,
               date_os     = date_enroll + day_os,
               date_event  = date_enroll + day_event)
}
