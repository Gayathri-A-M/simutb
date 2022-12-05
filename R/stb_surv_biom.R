## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating joint PFS and a biomarker.
##
## All functions start with "stb_surv_biom"
##
##
##
## -----------------------------------------------------------------------------

#' Theoretical derivation
#'
#' Theoretical derivation for selecting dose arm based on EFS with EFS being the
#' primary endpoint
#'
#' @export
#'
stb_surv_biom_theory <- function(drop_info_frac = 0.3,
                                 info_frac      = c(0.4, 1),
                                 n_arm          = 3,
                                 btype_primary  = "asOF",
                                 alpha          = 0.025,
                                 n_large        = 1000000,
                                 seed           = NULL,
                                 ...) {

    f_rej <- function(zs) {
        zs_other <- zs[1:n_other]
        zs_1     <- zs[n_other]
        zs_int   <- zs[-(1:n_other)]
        rej_1    <- all(zs_1 >= zs_other)
        rej_2    <- any(zs_int > boundary)
        rej_3    <- rej_1 & rej_2

        c(rej_1, rej_2, rej_3)
    }

    f_int <- function(zs) {
        rej <- f_rej(zs)
        rst <- rej[3] * dmvnorm(zs,
                                mean  = rep(0, length(zs)),
                                sigma = mat_sig)
        rst
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## get variance matrix
    n_int   <- length(info_frac)
    n_other <- n_arm - 1
    mat_sig <- matrix(NA, n_other + n_int, n_other + n_int)
    mat_sig[1:n_other, 1:n_other] <- 1 / 2
    for (k in 1:n_other) {
        for (j in 1:n_int) {
            if (k == n_other) {
                fac <- 1
            } else {
                fac <- 2
            }
            mat_sig[k, j + n_other] <- mat_sig[j + n_other, k] <-
                sqrt(drop_info_frac / info_frac[j]) / fac
        }
    }

    for (i in 1:n_int) {
        for (j in i:n_int) {
            mat_sig[i + n_other, j + n_other] <-
                mat_sig[j + n_other, i + n_other] <-
                sqrt(info_frac[i] / info_frac[j])
        }
    }

    diag(mat_sig) <- 1

    ## boundary
    boundary <-
        getDesignGroupSequential(sided            = 1,
                                 alpha            = alpha,
                                 informationRates = info_frac,
                                 typeOfDesign     = btype_primary)
    boundary <- boundary$criticalValues

    ## integration
    smps <- rmvnorm(n_large,
                    mean  = rep(0, nrow(mat_sig)),
                    sigma = mat_sig)

    rej <- apply(smps, 1, f_rej)
    rej <- t(rej)

    ## result
    rates <- c(margin      = mean(rej[, 2]),
               sel_arm     = mean(rej[, 1]),
               conditional = sum(rej[, 3]) / sum(rej[, 1]))

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    list(smps  = smps,
         rej   = rej,
         sigma = mat_sig,
         bound = boundary,
         rates = rates)
}


#' Theoretical derivation II
#'
#' Theoretical derivation for selecting dose arm based on pCR
#'
#' @param n sample size each arm at arm-selection look
#' @param p_biom true pcR rate
#' @param rho correlation between pCR statistic and log-rank test statistic
#' @param alpha alpha level at the primary analysis
#' @param margin margin to select arm 2
#'
#' @export
#'
stb_surv_biom_theory_pcr <- function(n_1, p_biom_1, rho_1,
                                     n_2        = n_1,
                                     p_biom_2   = p_biom_1,
                                     rho_2      = rho_1,
                                     margin     = 0.1,
                                     alpha      = 0.025,
                                     n_large    = 1000000,
                                     seed       = NULL,
                                     ...) {

    f_rej <- function(zs) {

        w <- zs[1:2]
        z <- zs[3:4]

        if (w[2] - w[1] > margin) {
            sel_arm <- 2
        } else {
            sel_arm <- 1
        }

        rej_arm          <- c(0, 0)
        rej_arm[sel_arm] <- z[sel_arm] > thresh

        c(sel_arm, rej_arm)
    }


    if (!is.null(seed))
        old_seed <- set.seed(seed)

    thresh  <- - qnorm(alpha)

    ## arm 1
    var_w   <- p_biom_1 * (1 - p_biom_1) / n_1
    rho_v   <- rho_1^2 * var_w
    w1      <- rnorm(n_large, p_biom_1, sqrt(var_w))
    z1      <- rnorm(n_large,
                     rho_1 * (w1 - p_biom_1),
                     sqrt(1 - rho_v))

    ## arm 2
    var_w   <- p_biom_2 * (1 - p_biom_2) / n_2
    rho_v   <- rho_2^2 * var_w
    w2      <- rnorm(n_large, p_biom_2, sqrt(var_w))
    z2      <- rnorm(n_large,
                     rho_2 * (w2 - p_biom_2),
                     sqrt(1 - rho_v))
    smps    <- cbind(w1, w2, z1, z2)

    ## rejection
    rej <- apply(smps, 1, f_rej)
    rej <- t(rej)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    sel_arm_1  <- mean(1 == rej[, 1])
    rej_arm_1  <- mean(rej[, 2])
    rej_arm_2  <- mean(rej[, 3])
    rej_any    <- rej_arm_1 + rej_arm_2
    rej_con_a1 <- mean(rej[which(1 == rej[, 1]), 2])
    rej_con_a2 <- mean(rej[which(2 == rej[, 1]), 3])

    c(p_biom_1   = p_biom_1,
      n_1        = n_1,
      rho_1      = rho_1,
      p_biom_2   = p_biom_2,
      n_2        = n_2,
      rho_2      = rho_2,
      sel_a1     = sel_arm_1,
      sel_a2     = 1 - sel_arm_1,
      rej_a1     = rej_arm_1,
      rej_a2     = rej_arm_2,
      rej_any    = rej_any,
      rej_con_a1 = rej_con_a1,
      rej_con_a2 = rej_con_a2)
}

#' Simulate an arm
#'
#' @param p_biomarker biomarker rates
#' @param median_mth median survival for responders and non-responsders
#' @export
#'
stb_surv_biom_arm_simu <- function(n,
                                   p_biomarker,
                                   median_mth,
                                   enroll_dur_mth,
                                   annual_drop = 0,
                                   date_bos    = as.Date("2022-1-1"),
                                   seed        = NULL,
                                   ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    smp_biom   <- rmultinom(1, n, p_biomarker)
    dta_enroll <- stb_tl_simu_enroll(n,
                                     enroll_dur_mth,
                                     date_bos    = date_bos,
                                     ...)
    dta_censor <- stb_tl_rexp(n,
                              median_mth  = NULL,
                              hazard      = NULL,
                              annual_drop = annual_drop,
                              ...)

    dta_biom  <- NULL
    dta_event <- NULL
    for (j in seq_len(length(smp_biom))) {
        cur_biom  <- rep(j - 1, smp_biom[j])
        cur_event <- stb_tl_rexp(smp_biom[j],
                                 median_mth  = median_mth[j],
                                 hazard      = NULL,
                                 annual_drop = NULL,
                                 ...)

        dta_biom  <- c(dta_biom,  cur_biom)
        dta_event <- c(dta_event, cur_event)
    }

    days <- stb_tl_pfs_os(dta_event,
                          dta_event,
                          dta_censor)

    rst <- cbind(dta_enroll, days, biom = dta_biom) %>%
        data.frame() %>%
        select(-day_prog, -day_dth, -day_event) %>%
        mutate(date_censor = date_enroll + day_censor,
               date_pfs    = date_enroll + day_pfs,
               date_os     = date_enroll + day_os)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Simulate a trial
#'
#'
#' @export
#'
stb_surv_biom_trial_simu <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n_by_arm    <- rmultinom(1,
    ##                          lst_design$sample_size,
    ##                          lst_design$ratio_by_arm)
    n_by_arm           <- lst_design$sample_size
    n_by_arm           <- n_by_arm * lst_design$ratio_by_arm
    n_by_arm           <- n_by_arm / sum(lst_design$ratio_by_arm)
    n_by_arm           <- floor(n_by_arm)

    p_biom_by_arm      <- lst_design$p_biom_by_arm
    median_surv_by_arm <- lst_design$median_surv_by_arm

    rst <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_arm <- stb_surv_biom_arm_simu(
            n              = n_by_arm[i],
            p_biomarker    = p_biom_by_arm[i, ],
            median_mth     = median_surv_by_arm[i, ],
            enroll_dur_mth = lst_design$enroll_dur_mth,
            annual_drop    = lst_design$annual_drop,
            ...)

        cur_arm$arm <- as.character(i - 1)
        rst         <- rbind(rst, cur_arm)
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Plot a simulated data
#'
#'
#' @export
#'
stb_surv_biom_trial_plot <- function(data, facet_by = "arm") {
    fit      <- survfit(Surv(day_pfs, status_pfs) ~ arm + biom, data = data)
    plt_pfs  <- ggsurvplot_facet(fit,
                                 data = data,
                                 facet.by = c(facet_by))

    plt_pfs
}

#' Arm Selection Rule I
#'
#' Arm selection by interim biomarker. Select the treatment arm with the largest
#' biomarker mean.
#'
#' @export
#'
stb_surv_biom_arm_sel_rule_1 <- function(data, ...) {
    dta_biom <- data %>%
        filter(arm > 0) %>%
        group_by(arm) %>%
        summarize(n    = n(),
                  biom = mean(biom)) %>%
        arrange(biom)

    sel_arm <- as.numeric(dta_biom[nrow(dta_biom),
                                   'arm'])

    list(sel_arm = sel_arm,
         biom    = dta_biom)
}

#' Arm Selection Rule II
#'
#' Arm selection by pfs. Select the treatment arm with the best
#' hazard rate.
#'
#' @export
#'
stb_surv_biom_arm_sel_rule_2 <- function(data, fml_surv) {
    test_pfs <- stb_tl_surv_logrank(data = data,
                                    fml  = fml_surv)

    zscore  <- which(grepl("zscore", names(test_pfs)))
    zscore  <- test_pfs[zscore]
    sel_arm <- which.min(zscore)

    list(sel_arm = sel_arm,
         biom    = data.frame(arm  = seq_len(length(zscore)),
                              n    = 0, ## place holder
                              biom = zscore))
}

#' Arm Selection Rule I
#'
#' Arm selection by interim biomarker. Select arm 2 only when its biom rate is
#' xx% better than arm 1
#'
#' @export
#'
stb_surv_biom_arm_sel_rule_3 <- function(data, margin = 0.1, ...) {
    dta_biom <- data %>%
        filter(arm > 0) %>%
        group_by(arm) %>%
        summarize(n    = n(),
                  biom = mean(biom)) %>%
        arrange(arm)

    rate_arm_1 <- dta_biom[1, "biom"]
    rate_arm_2 <- dta_biom[2, "biom"]

    if (rate_arm_2 - rate_arm_1 > margin) {
        sel_arm <- 2
    } else {
        sel_arm <- 1
    }

    list(sel_arm = sel_arm,
         biom    = dta_biom)
}


#' Simulate a trial with interim
#'
#'
#'
#' @export
#'
stb_surv_biom_trial_interim <- function(lst_design, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    data_full <- stb_surv_biom_trial_simu(lst_design)
    if (0) {
        ## check
        data_full %>%
            group_by(arm) %>%
            summarize(biom = mean(biom))

        data_full %>%
            group_by(arm, biom) %>%
            summarize(n = n(),
                      pfs = median(day_pfs) / 30.1)
    }

    ## interim look data for arm selection
    data_drop <- stb_tl_interim_data(data_full,
                                     total     = lst_design$drop_target,
                                     info_frac = lst_design$drop_info_frac,
                                     event     = lst_design$drop_event)

    ## arm selection
    f_sel        <- lst_design$f_arm_sel
    lst_arm_sel  <- f_sel(data_drop,
                          lst_design$fml_surv,
                          margin = lst_design$margin_r3)

    ## data selected
    data_sel <- data_full %>% filter(arm %in% c(0, lst_arm_sel$sel_arm))

    ## interim analysis
    info_frac <- lst_design$info_frac
    rst       <- NULL
    for (a in 1:2) {
        data_sel <- data_full %>%
            filter(arm %in% c(0, a))

        for (i in seq_len(length(lst_design$info_frac))) {
            data_interim <- stb_tl_interim_data(
                data_sel,
                total     = lst_design$target_events,
                info_frac = info_frac[i],
                event     = "pfs")

            test_pfs <- stb_tl_surv_logrank(data = data_interim,
                                            fml  = lst_design$fml_surv)

            pval     <- unname(test_pfs["pvalue_oneside"])
            bound    <- lst_design$bound[i]
            rst      <- rbind(rst,
                              c(sel_arm    = unname(lst_arm_sel$sel_arm),
                                arm        = a,
                                inx        = i,
                                info_frac  = info_frac[i],
                                bound      = bound,
                                pval       = pval,
                                rej        = pval <= bound,
                                zscore     = unname(test_pfs["zscore"]),
                                nevent     = unname(test_pfs["nevent"]),
                                hr         = unname(test_pfs["hr"])))

        }
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rownames(rst) <- NULL
    list(rst     = data.frame(rst),
         biom    = lst_arm_sel$biom)
}


#' Summarize simulation results
#'
#' @export
#'
stb_surv_biom_summary <- function(lst_rst) {

    n_reps      <- length(lst_rst)
    rst_biom    <- list()
    rst_interim <- list()
    for (i in seq_len(length(lst_rst))) {
        rst_biom[[i]]    <- lst_rst[[i]]$biom %>% mutate(rep = i)
        rst_interim[[i]] <- lst_rst[[i]]$rst  %>% mutate(rep = i)
    }

    rst_biom    <- rbindlist(rst_biom)
    rst_interim <- rbindlist(rst_interim)

    ## arm selection
    rst_sel_arm <- rst_interim %>%
        distinct(sel_arm, rep) %>%
        group_by(sel_arm) %>%
        summarize(n    = n(),
                  rate = n / n_reps)

    ## biomarker at interim
    rst_biom_summary <- rst_biom %>%
        group_by(arm) %>%
        summarize(n    = mean(n),
                  mean = mean(biom))

    ## estimate
    rst_hr <- rst_interim %>%
        group_by(arm) %>%
        summarize(hr = mean(hr))

    ## rejection
    rej <- rst_interim %>%
        filter(sel_arm == arm &
               rej == 1) %>%
        group_by(rep, sel_arm) %>%
        arrange(inx, .by_group = TRUE) %>%
        slice(n = 1) %>%
        ungroup() %>%
        group_by(inx, info_frac, sel_arm) %>%
        summarize(N_Rej = n()) %>%
        mutate(Rej = N_Rej / n_reps) %>%
        group_by(sel_arm) %>%
        mutate(CumuRej = cumsum(Rej))

    rej_any <- rej %>%
        group_by(inx, info_frac) %>%
        summarize(N_Rej   = sum(N_Rej),
                  Rej     = sum(Rej),
                  CumuRej = sum(CumuRej))

    ## return
    rst <- list(sel_arm     = rst_sel_arm,
                biom        = rst_biom_summary,
                rej_arm     = rej,
                rej_any     = rej_any,
                hr          = rst_hr,
                raw_interim = rst_interim,
                raw_biom    = rst_biom)

    rst
}

#' Key simulation results
#'
#' @export
#'
stb_surv_biom_key <- function(lst_design, rst_summary, seed) {

    fp <- function(vname) {
        paste(lst_design[[vname]],
              collapse = ",")
    }

    frej <- function(dta, a = NULL) {
        dta <- dta %>%
            arrange(info_frac)

        if (!is.null(a)) {
            dta <- dta %>%
                filter(sel_arm == a)
        }

        unname(as.numeric(dta[nrow(dta), "CumuRej"]))
    }

    rst_key <- data.frame(
        info_frac      = fp("info_frac"),
        hr             = fp("hr"),
        ctl_median     = lst_design$ctl_median_surv[2],
        sample_size    = lst_design$sample_size,
        target_events  = lst_design$target_events,
        drop_event     = lst_design$drop_event,
        drop_info_frac = lst_design$drop_info_frac,
        btype_primary  = lst_design$btype_primary,
        sel_arm_2      = unname(as.numeric(rst_summary$sel_arm[2, "rate"])),
        rej_any        = frej(rst_summary$rej_any),
        rej_arm_1      = frej(rst_summary$rej_arm, 1),
        rej_arm_2      = frej(rst_summary$rej_arm, 2),
        seed           = seed)

    rst_key
}
