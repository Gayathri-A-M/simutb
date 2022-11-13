## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating stratified survival analysis.
## All functions start with "stb_surv_stra"
##
##
##
## -----------------------------------------------------------------------------



#' Simulate an arm
#'
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
stb_surv_stra_stra_simu <- function(n,
                                     median_mth,
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

    dta_event  <- stb_tl_rexp(n,
                              median_mth  = median_mth,
                              hazard      = NULL,
                              annual_drop = NULL,
                              mth_to_days = mth_to_days)

    dta_censor <- stb_tl_rexp(n,
                              median_mth  = NULL,
                              hazard      = NULL,
                              annual_drop = annual_drop,
                              mth_to_days = mth_to_days)

    days <- stb_tl_pfs_os(floor(dta_event),
                          floor(dta_event),
                          floor(dta_censor))

    cbind(dta_enroll, days) %>%
        data.frame() %>%
        select(-day_prog, -day_dth, -day_event) %>%
        mutate(date_censor = date_enroll + day_censor,
               date_pfs    = date_enroll + day_pfs,
               date_os     = date_enroll + day_os)
}

#' Simulate a trial
#'
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
stb_surv_strat_trial_simu <- function(lst_design,
                                      seed = NULL, ...) {
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    n_stra      <- length(lst_design$stra_freq)
    median_surv <- rbind(lst_design$ctl_median,
                         lst_design$ctl_median / lst_design$hr)
    stra_size   <- rmultinom(1,
                             lst_design$sample_size,
                             lst_design$stra_freq)

    arm_ratio   <- c(lst_design$ctl_ratio,
                     1 - lst_design$ctl_ratio)

    rst <- NULL
    for (j in 1:n_stra) {
        cur_stra_size  <- stra_size[j]
        for (arm in 1:2) {
            cur_arm <- floor(cur_stra_size * arm_ratio[arm])
            cur_rst <- stb_surv_stra_stra_simu(
                n              = cur_arm,
                median_mth     = median_surv[arm, j],
                enroll_dur_mth = lst_design$enroll_dur_mth,
                annual_drop    = lst_design$annual_drop,
                ...)

            cur_rst$arm     <- arm - 1
            cur_rst$stratum <- j
            rst             <- rbind(rst, cur_rst)
        }
    }


    ## keep the first events
    rst <- stb_tl_interim_data_2arm(rst, lst_design$target, 1)

    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}
