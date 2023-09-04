## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for simulating enrollment time
##
##  DATE:
##      APRIL, 2023
## -----------------------------------------------------------------------------


#' Simulate Enrollment Time
#'
#' An archived version of enrollment simulation (archived function)
#'
#' @export
#'
stb_tl_simu_enroll_arc <- function(ntot,
                                   enroll_dur_mth,
                                   min_fu      = NULL,
                                   date_bos    = NULL,
                                   mth_to_days = 30.4,
                                   ...) {

    rand_enroll <- runif(ntot, 0, enroll_dur_mth)
    day_enroll  <- ceiling(rand_enroll * mth_to_days)

    rst <- data.frame(sid        = seq_len(ntot),
                      day_enroll = day_enroll)

    ## set up dates in addition to days
    if (!is.null(date_bos)) {
        rst$date_bos    <- date_bos
        rst$date_enroll <- date_bos + day_enroll
    }

    ## set up end of study time by minimum follow up
    if (!is.null(min_fu)) {
        day_eos      <- max(rand_enroll) + min_fu
        day_eos      <- day_eos * mth_to_days
        day_eos      <- floor(day_eos)
        rst$day_eos  <- day_eos - day_enroll
        rst$date_eos <- date_bos + day_eos
    }

    ## return
    rst
}


#' Simulate Enrollment Time
#'
#' Simulate enrollment time by total time
#'
#' @export
#'
stb_tl_simu_enroll_by_dur <- function(n_pt, n_pt_tot, pt_dur_mth) {
    data.frame(mth_enroll = runif(n_pt,
                                  0,
                                  pt_dur_mth))
}

#' Simulate Enrollment Time
#'
#' Simulate enrollment time by rate
#'
#' @export
#'
stb_tl_simu_enroll_by_rate <- function(n_pt, n_pt_tot, pt_per_mth) {
    pt_dur_mth <- n_pt_tot / pt_per_mth
    data.frame(mth_enroll = runif(n_pt,
                                  0,
                                  pt_dur_mth))
}

#' Simulate Enrollment Time
#'
#' Simulate enrollment time by center
#'
#' @export
#'
stb_tl_simu_enroll_by_center <- function(n_pt,
                                         n_pt_tot,
                                         n_center,
                                         pt_per_center_per_mth,
                                         center_per_mth) {

    center_dur_mth <- n_center / center_per_mth
    pt_dur_mth     <- n_pt_tot / pt_per_center_per_mth

    en_center      <- runif(n_center, 0, center_dur_mth)
    en_center      <- sort(en_center)

    rst <- NULL
    for (i in 1:n_center) {
        cur_center <- en_center[i]
        cur_en_pt  <- cur_center + runif(n_pt, 0, pt_dur_mth)
        rst        <- rbind(rst, cbind(center            = i,
                                       mth_enroll_center = cur_center,
                                       mth_enroll        = cur_en_pt))
    }

    data.frame(rst) %>%
        arrange(mth_enroll) %>%
        slice_head(n = n_pt)
}

#' Simulate Enrollment Time
#'
#' @param mth_fix_fu fixed follow-up months for all patients
#' @param mth_min_fu minimum follow-up months for all patients
#'
#' @export
#'
stb_tl_simu_enroll <- function(n_pt_arm    = 100,
                               n_pt_tot    = n_pt_arm * 2,
                               par_enroll  = list(type       = "by_duration",
                                                  pt_dur_mth = 24),
                               mth_min_fu  = NULL,
                               mth_fix_fu  = NULL,
                               date_bos    = "2022-01-01",
                               mth_to_days = 30.4,
                               ...) {

    type                <- par_enroll$type
    par_enroll$type     <- NULL
    par_enroll$n_pt     <- n_pt_arm
    par_enroll$n_pt_tot <- n_pt_tot
    f_enroll            <- switch(type,
                                  by_duration = stb_tl_simu_enroll_by_dur,
                                  by_rate     = stb_tl_simu_enroll_by_rate,
                                  by_center   = stb_tl_simu_enroll_by_center)

    rst <- do.call(f_enroll, par_enroll) %>%
        mutate(day_enroll = mth_to_days * mth_enroll) %>%
        arrange(day_enroll) %>%
        mutate(sid = row_number())

    ## set up dates in addition to days
    if (!is.null(date_bos)) {
        rst$date_bos    <- as.Date(date_bos)
        rst$date_enroll <- as.Date(date_bos) + rst$day_enroll
    }

    ## set up end of study time by minimum or fixed follow up
    day_chk <- NULL
    if (!is.null(mth_min_fu)) {
        day_chk <- max(rst$day_enroll)
        min_fu  <- mth_min_fu
    } else if (!is.null(mth_fix_fu)) {
        day_chk <- rst$day_enroll
        min_fu  <- mth_fix_fu
    }

    if (!is.null(day_chk)) {
        day_eos      <- day_chk + min_fu * mth_to_days
        day_eos      <- floor(day_eos)
        rst$day_eos  <- day_eos - rst$day_enroll
        rst$date_eos <- rst$date_bos + day_eos
    }

    ## return
    rst
}

#' Simulate Enrollment Time
#'
#' Simulate Enrollment Time for all Arms
#'
#' @export
#'
stb_tl_simu_enroll_arms <- function(n_by_arm, ..., seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    n_arm <- length(n_by_arm)
    rst   <- NULL
    for (i in seq_len(n_arm)) {
        cur_arm     <- stb_tl_simu_enroll(n_by_arm[i],
                                          n_pt_tot = sum(n_by_arm),
                                          ...)
        cur_arm$arm <- i - 1
        rst         <- rbind(rst, cur_arm)
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}
