## -----------------------------------------------------------------------------
##
##     DESCRIPTION:
##
##          This file contains R functions for simulating bivariate survival
##          time, and borrow historical log HR information.
##
##     DATE:
##          DEC, 2023
##
##     BY:
##          STAT INNOVATION
##
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
surv_bor_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Time to event Bayesian analysis with historical information \n")
    cat("    borrowing. \n\n")

    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 500)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(1/2, 1/2))\n")
    cat("    ctl_median_mth: control group median survival time in months\n")
    cat("                    for each endpoint. default c(12, 12) \n")
    cat("    hr:             hazard ratio (treatment vs. control) for\n")
    cat("                    each endpoint. default c(0.7, 0.7)\n")
    cat("    target_primary: target number of event for the first endpoint\n")
    cat("                    (default 20)\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    date_bos:       begin of study date\n")
    cat("                    (default 2022-1-1) \n")
    cat("    par_analysis:   list of analysis parameters\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}

#' Set design paramter
#'
#'
surv_bor_default_para <- function(x) {
    lst_default <- list(
        sample_size     = 500,
        ratio_by_arm    = c(1, 1) / 2,
        ctl_median_mth  = 12,
        hr              = 0.7,
        target_primary  = 20,
        annual_drop     = 0.000001,
        mth_fix_fu      = 600,

        par_enroll = list(
            type       = "by_duration",
            pt_dur_mth = 12),

        par_analysis = list(
            prior   = c(loghr_mean = -0.7,
                        loghr_sd   = 1),
            loghr_threshold = 0,
            a               = 1)
        )
}

#' Simulate a trial
#'
#'
#' Generate data
#'
#'
#'
surv_bor_gen_data <- function(lst_design,
                              seed = NULL,
                              ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n by arm
    n_by_arm <- tl_draw_arm_size(sample_size  = lst_design$sample_size,
                                 ratio_by_arm = lst_design$ratio_by_arm)

    ## enrollment
    dat_enroll <- do.call(stb_tl_simu_enroll_arms,
                          list(n_by_arm   = c(n_by_arm),
                               mth_fix_fu = lst_design$mth_fix_fu,
                               par_enroll = lst_design$par_enroll))

    ## outcome
    rst <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_censor <- stb_tl_rexp(n_by_arm[i],
                                  median_mth  = NULL,
                                  hazard      = NULL,
                                  annual_drop = lst_design$annual_drop)
        if (1 == i) {
            median_mth <- lst_design$ctl_median_mth
        } else {
            median_mth <- lst_design$ctl_median_mth * lst_design$hr
        }

        cur_event <- stb_tl_rexp(n_by_arm[i],
                                 median_mth  = median_mth,
                                 hazard      = NULL,
                                 annual_drop = NULL)

        cur_rst <- data.frame(arm        = i - 1,
                              sid        = 1:n_by_arm[i],
                              day_censor = cur_censor,
                              day_event  = cur_event)

        rst <- rbind(rst, cur_rst)
    }

    ## merge data
    rst <- dat_enroll %>%
        left_join(rst, by = c("arm", "sid")) %>%
        mutate(date_censor = date_enroll + day_censor,
               date_event  = date_enroll + day_event) %>%
        rowwise() %>%
        mutate(day_obs    = min(day_event, day_censor, day_eos),
               status_obs = if_else(day_event == day_obs, 1, 0),
               date_obs   = date_enroll + day_obs) %>%
        ungroup()


    ## make target events
    if (!is.null(lst_design$target_primary)) {
        dat_target <- rst %>%
            filter(1 == status_obs) %>%
            arrange(date_obs)

        if (nrow(dat_target) > lst_design$target_primary) {
            date_target <- dat_target$date_obs[lst_design$target_primary]
            rst <- rst %>%
                mutate(date_target = date_target) %>%
                mutate(date_obs_old = date_obs,
                       status_obs   = if_else(date_obs_old <= date_target,
                                              status_obs,
                                              0),
                       date_obs = if_else(date_obs_old <= date_target,
                                          date_obs_old,
                                          date_target),
                       day_obs  = date_obs - date_enroll)
        }
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Borrow from historical information
#'
#'
surv_bor_ana_post <- function(loghr_mean, loghr_sd,
                             prior_loghr_mean, prior_loghr_sd,
                             a = 0) {

    tmp     <- a / prior_loghr_sd^2 + 1 / loghr_sd^2
    mu_post <- a * prior_loghr_mean / prior_loghr_sd^2
    mu_post <- mu_post + loghr_mean / loghr_sd^2

    mu_post <- mu_post / tmp
    sd_post <- sqrt(1 / tmp)

    c(mu_post = unname(mu_post),
      sd_post = unname(sd_post))
}


#' Analyze a trial
#'
#'
#'
surv_bor_ana_data <- function(data_full, par_analysis, n_post = 3000, ...) {
    cox_fit <- stb_tl_surv_logrank(
        data_full,
        fml_survdiff = "Surv(day_obs, status_obs) ~ arm")

    post_dist <- surv_bor_ana_post(
        cox_fit["log_hr"],
        cox_fit["se_log_hr"],
        par_analysis$prior["loghr_mean"],
        par_analysis$prior["loghr_sd"],
        a = par_analysis$a)

    post_smps <- rnorm(n_post,
                       post_dist["mu_post"],
                       post_dist["sd_post"])


    pbigger <- mean(post_smps > par_analysis$loghr_threshold)

    ## return
    data.frame(
        loghr_mean = cox_fit["log_hr"],
        loghr_sd   = cox_fit["se_log_hr"],
        pval_cox   = cox_fit["pval_cox"],
        nevent     = cox_fit["nevent"],
        post_loghr_mean = post_dist["mu_post"],
        post_loghr_sd   = post_dist["sd_post"],
        pbigger         = pbigger)
}
