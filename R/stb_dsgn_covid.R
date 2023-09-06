## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##
##      This file contains R functions for multi-stage multi-arm study
##      with survival endpoint
##
##  DATE:
##
##      MAY, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
covid_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Multi-Stage Multi-Arm Trial with Survival Outcome\n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(1/3, 1/3, 1/3))\n")
    cat("    six_mth_cumu:   six month cumulative incidence rate for\n")
    cat("                    each arm (default c(0.021, 0.005, 0.005))\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    min_fix_fu:     fixed follow-up months (default 6) \n")
    cat("    target_primary: target primary event (default 54) \n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    power:          pairwise comparison power (default 0.9)\n")
    cat("    prior_weight:   weight on the prior knowledge \n")
    cat("                    (default 0.05)\n")
    cat("    par_interim:    list of parameters for creating interim \n")
    cat("                    and final analysis dataset\n")
    cat("    par_analysis:   list of analysis parameters\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}

#' Default design parameter
#'
#'
internal_covid_dpara <- function(lst) {
    lst$median_mth     <- NULL
    lst$target_primary <- NULL
    lst$six_mth_cumu   <- c(0.021, 0.005, 0.005)
    lst$prior_weight   <- 0.05
    lst$par_interim    <- list(target_primary = 54,
                               ana_fraction   = c(0.5, 1),
                               target_arms    = NULL)
    lst$power          <- 0.9

    lst$par_analysis   <- list(futility_cond_pow   = 0.4,
                               futility_pred_pow   = 0.4,
                               prior_alpha         = 0.01,
                               prior_beta          = 0.01,
                               futility_ctl_fu     = 12,
                               futility_ctl_rate   = 0.02,
                               futility_ctl_rate_c = 0.5)

    lst
}


#' Generate analysis data set
#'
#' Generate analysis dataset for minimum fu days
#'
covid_surv_ana_logrank <- function(data_ana,
                                   fml = "Surv(day_obs, status_obs) ~ arm",
                                   mth_fix_fu = NULL) {

    day_analysis <- data_ana %>%
        mutate(day = date_interim - date_bos) %>%
        select(day)


    n_arm <- max(data_ana$arm)
    rst   <- NULL
    for (i in 1:n_arm) {

        cur_data <- data_ana %>%
            filter(arm %in% c(0, i)) %>%
            mutate(arm = if_else(arm > 0, 1, 0))

        cur_rst   <- stb_tl_surv_logrank(cur_data,
                                         fml_survdiff = fml)

        n_enroll  <- cur_data %>%
            group_by(arm) %>%
            summarize(n = n(),
                      y = sum(day_obs))

        n_event   <- cur_data %>%
            filter(1 == status_obs) %>%
            group_by(arm) %>%
            summarize(n = n())

        n_fu <- c(NA, NA)
        if (!is.null(mth_fix_fu)) {
            n_fu <- cur_data %>%
                mutate(day_on_study = date_interim - date_enroll) %>%
                filter(day_on_study >= mth_fix_fu * 30.5)

            if (0 == nrow(n_fu)) {
                n_fu <- c(0, 0)
            } else {
                n_fu <- n_fu %>%
                    group_by(arm) %>%
                    summarize(n = n())

                n_fu <- n_fu$n
            }
        }

        rst <- rbind(rst,
                     c(arm = i,
                       n_enroll_ctl = n_enroll$n[1],
                       n_enroll_trt = n_enroll$n[2],
                       n_event_ctl  = n_event$n[1],
                       n_event_trt  = n_event$n[2],
                       n_toty_ctl   = n_enroll$y[1],
                       n_toty_trt   = n_enroll$y[2],
                       n_event      = sum(n_event$n),
                       n_min_fu_ctl = n_fu[1],
                       n_min_fu_trt = n_fu[2],
                       n_min_fu     = sum(n_fu),
                       n_enroll     = sum(n_enroll$n),
                       day_analysis = day_analysis$day[1],
                       cur_rst["pval_oneside"],
                       cur_rst["hr"],
                       cur_rst["se_log_hr"],
                       cur_rst["zscore"]))
    }

    data.frame(rst)
}


#' Analyze interim and final data
#'
#'
covid_analysis <- function(dta_ana, info_frac, lst_design) {

    cur_rst <- covid_surv_ana_logrank(
        dta_ana,
        mth_fix_fu = lst_design$mth_fix_fu)

    cur_rst         <- data.frame(cur_rst)
    cur_rst$interim <- info_frac

    ## conditional power
    cur_rst$cond_pow <- stb_tl_gsd_condpower(
        cur_rst$zscore,
        info_frac    = as.numeric(info_frac),
        alpha        = lst_design$alpha,
        power        = lst_design$power,
        use_observed = FALSE)

    cur_rst$cond_pow_hat <- stb_tl_gsd_condpower(
        cur_rst$zscore,
        info_frac    = as.numeric(info_frac),
        alpha        = lst_design$alpha,
        power        = lst_design$power,
        use_observed = TRUE)

    ## predictive power
    cur_rst$pred_pow <- stb_tl_gsd_predpower(
        cur_rst$zscore,
        info_frac    = as.numeric(info_frac),
        alpha        = lst_design$alpha,
        power        = lst_design$power,
        prior_weight = lst_design$prior_weight)


    ## return
    cur_rst
}

#' Summarize the simulation results
#'
#'
covid_simu_summary_ind <- function(dat_rst, par_analysis, nsmps = 3000) {

    futility_cond_pow    <- par_analysis$futility_cond_pow
    futility_pred_pow    <- par_analysis$futility_pred_pow
    prior_alpha          <- par_analysis$prior_alpha
    prior_beta           <- par_analysis$prior_beta
    futility_ctl_fu      <- par_analysis$futility_ctl_fu * 30.5
    futility_ctl_rate    <- par_analysis$futility_ctl_rate
    futility_ctl_rate_c  <- par_analysis$futility_ctl_rate_c

    dat_rst <- dat_rst %>%
        mutate(interim = as.numeric(interim))

    n_arm      <- max(dat_rst$arm)
    info_fracs <- sort(unique(dat_rst$interim))
    rst        <- NULL
    for (carm in seq_len(n_arm)) {
        for (cinfo in seq_len(length(info_fracs))) {
            cur_rst <- dat_rst %>%
                dplyr::filter(arm     == carm &
                              interim == info_fracs[cinfo])

            post_alpha      <- prior_alpha - 1 + cur_rst$n_event_ctl[1]
            post_beta       <- prior_beta  + cur_rst$n_toty_ctl[1]
            post_smps       <- rgamma(nsmps, post_alpha, post_beta)
            post_rates      <- 1 - exp(-post_smps * futility_ctl_fu)

            futile_cond     <- cur_rst$cond_pow[1]     < futility_cond_pow
            futile_cond_hat <- cur_rst$cond_pow_hat[1] < futility_cond_pow
            futile_pred     <- cur_rst$pred_pow[1]     < futility_pred_pow
            post_ctl_rate   <- mean(post_rates)
            prob_ctl_rate   <- mean(post_rates > futility_ctl_rate)
            futile_ctl_rate <- prob_ctl_rate < futility_ctl_rate_c

            rst <- rbind(rst,
                         data.frame(arm             = carm,
                                    interim         = cinfo,
                                    info_frac       = info_fracs[cinfo],
                                    hr              = cur_rst$hr[1],
                                    day_analysis    = cur_rst$day_analysis[1],
                                    post_ctl_rate   = post_ctl_rate,
                                    prob_ctl_rate   = prob_ctl_rate,
                                    futile_cond     = futile_cond,
                                    futile_cond_hat = futile_cond_hat,
                                    futile_pred     = futile_pred,
                                    futile_ctl_rate = futile_ctl_rate))
        }
    }

    rst
}


#' Summarize the simulation results
#'
#'
#' @export
#'
covid_simu_summary <- function(lst, par_analysis) {
    rst <- NULL
    for (i in seq_len(length(lst))) {
        cur_rst <- covid_simu_summary_ind(lst[[i]][[1]],
                                          par_analysis)
        cur_rst$rep <- i
        rst         <- rbind(rst, cur_rst)
    }

    list(rst)
}
