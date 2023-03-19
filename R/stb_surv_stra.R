## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating stratified survival analysis.
## All functions start with "stb_surv_stra"
##
##
##
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
stratsurv_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Stratified Survival Analysis\n\n")
    cat("Design Parameters:\n")
    cat("    target:         target number of events (default 100) \n")
    cat("    sample_size:    total sample size for treatment and control\n")
    cat("                    (default 500)\n")
    cat("    ctl_ratio:      randomization proportion for control \n")
    cat("                    (default 0.5)\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    enroll_dur_mth: enrollment months (default 18) \n")
    cat("    stra_freq:      frequency for each stratum \n")
    cat("                    (default c(0.5, 0.5)) \n")
    cat("    ctl_median:     control arm median survival months in \n")
    cat("                    each stratum (default c(15, 15)) \n")
    cat("    hr:             hazard ratio in each stratum\n")
    cat("                    (default c(0.7, 0.7)) \n")
    cat("    date_bos:       begin of study date\n")
    cat("                    (default 2022-1-1) \n")
}

#' Default design paramters
#'
#'
stratsurv_default_para <- function() {
    list(target          = 100,
         sample_size     = 500,
         ctl_ratio       = 0.5,
         annual_drop     = 0.000001,
         enroll_dur_mth  = 18,
         stra_freq       = c(0.5, 0.5),
         ctl_median      = c(15, 15),
         hr              = c(0.7, 0.7),
         mth_to_days     = 30.4,
         date_bos        = as.Date("2022-1-1"))
}

#' Simulate an arm
#'
#'
strasurv_simu_arm <- function(n,
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

#' Generate data for one trial
#'
#'
strasurv_gen_data <- function(lst_design,
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
            cur_rst <- strasurv_simu_arm(
                n              = cur_arm,
                median_mth     = median_surv[arm, j],
                enroll_dur_mth = lst_design$enroll_dur_mth,
                annual_drop    = lst_design$annual_drop,
                date_bos       = lst_design$date_bos,
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
