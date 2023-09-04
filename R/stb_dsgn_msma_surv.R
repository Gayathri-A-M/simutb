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
#' @export
#'
msma_surv_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Multi-Stage Multi-Arm Trial with Survival Outcome\n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(0.5, 0.5))\n")
    cat("    median_mth:     median survival in months for each arm\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    min_fix_fu:     fixed follow-up months (default 6) \n")
    cat("    target_primary: target primary event (default 54) \n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    par_interim:    list of parameters for creating interim \n")
    cat("                    and final analysis dataset\n")
    cat("    par_analysis:   list of analysis parameters\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}


#' Default design parameter
#'
#'
internal_msma_surv_dpara <- function() {
    list(sample_size        = 4800,
         ratio_by_arm       = c(1, 1, 1) / 3,
         median_mth         = c(10, 10, 10),
         annual_drop        = 0.03,
         par_enroll         = list(type       = "by_duration",
                                   pt_dur_mth = 16),
         mth_fix_fu         = 6,
         target_primary     = 54,
         par_interim        = cbind(ana_fraction = c(0.5, 1)),
         par_analysis       = list(),
         alpha              = 0.025
         )
}

#' Generate data
#'
#'
#'
msma_surv_gen_data <- function(lst_design,
                               seed = NULL,
                               set_target = FALSE,
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

    ## convert to median survival
    if (is.null(lst_design$median_mth) &
        !is.null(lst_design$six_mth_cumu))
        lst_design$median_mth <- stb_tl_cumu_median(lst_design$six_mth_cumu)

    ## outcome
    rst <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_censor <- stb_tl_rexp(n_by_arm[i],
                                  median_mth  = NULL,
                                  hazard      = NULL,
                                  annual_drop = lst_design$annual_drop)

        cur_event <- stb_tl_rexp(n_by_arm[i],
                                 median_mth  = lst_design$median_mth[i],
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


#' Generate analysis data set
#'
#' Generate analysis dataset for minimum fu days
#'
msma_surv_ana_logrank <- function(data_ana,
                                  fml = "Surv(day_obs, status_obs) ~ arm",
                                  mth_fix_fu = NULL) {

    n_arm <- max(data_ana$arm)
    rst   <- NULL
    for (i in 1:n_arm) {

        cur_data <- data_ana %>%
            filter(arm %in% c(0, i)) %>%
            mutate(arm = if_else(arm > 0, 1, 0))

        cur_rst   <- stb_tl_surv_logrank(cur_data,
                                         fml_survdiff = fml)

        rst <- rbind(rst,
                     c(arm = i,
                       cur_rst["pval_oneside"],
                       cur_rst["hr"],
                       cur_rst["zscore"]))
    }

    data.frame(rst)
}


#' Summarize simulation study results
#'
#'
msma_surv_simu_summary <- function(lst, ...) {
    rst <- NULL
    for (i in 1:length(lst)) {
        cur_rst     <- lst[[i]][[1]]
        cur_rst$rep <- i
        rst         <- rbind(rst, cur_rst)
    }

    list(rst)
}
