## -----------------------------------------------------------------------------
##
## This file contains R tool functions for survival analysis.
##
## All functions start with "stb_tl_".
##
##
##
## -----------------------------------------------------------------------------


#' Convert median survival or annual_drop to hazard
#'
#' @export
#'
stb_tl_hazard <- function(median_surv = NULL, annual_drop = NULL) {

    rst <- NULL
    if (!is.null(median_surv))
        rst <- c(rst,
                 - log(0.5) / median_surv)

    if (!is.null(annual_drop))
        rst <- c(rst,
                 - log(1 - annual_drop) / 12)

    rst
}


#' Simulate time to events in days
#'
#' @param ntot  total number of patients
#' @param hazard hazard for event
#' @param median_mth median survival in months
#' @param annual_drop annual drop rate
#'
#' @export
#'
stb_tl_rexp <- function(ntot,
                        median_mth   = NULL,
                        hazard       = NULL,
                        annual_drop  = NULL,
                        mth_to_days  = 30.4,
                        take_floor   = TRUE) {

    if (is.null(hazard)) {
        hazard <- stb_tl_hazard(median_surv = median_mth,
                                annual_drop = annual_drop)[1]
    }

    rand_event <- rexp(ntot, hazard) * mth_to_days

    if (take_floor)
        rand_event <- floor(rand_event)

    rand_event
}

#' Get PFS and OS
#'
#'
#'
#' @export
#'
stb_tl_pfs_os <- function(day_prog, day_dth, day_censor = NULL) {
    f_s <- function(prog, dth, censor) {
        day_pfs <- min(prog, dth, censor)
        day_os  <- min(dth, censor)

        status_pfs <- censor > day_pfs
        status_os  <- censor > day_os

        day_event <- min(prog, dth)

        c(prog,
          dth,
          censor,
          day_pfs,
          status_pfs,
          day_os,
          status_os,
          day_event)
    }

    if (is.null(day_censor))
        day_censor <- rep(Inf, length(day_prog))

    rst <- apply(cbind(day_prog, day_dth, day_censor),
                 1,
                 function(x) f_s(x[1], x[2], x[3]))

    rst <- t(rst)
    colnames(rst) <- c("day_prog",
                       "day_dth",
                       "day_censor",
                       "day_pfs",
                       "status_pfs",
                       "day_os",
                       "status_os",
                       "day_event")

    rst
}

#' Simulate Enrollment Time
#'
#' @inheritParams stb_simu_trial_single
#'
#' @param date_bos date begin of study
#' @param enroll_dur_mth enroll duration in months
#' @param min_fu minimum follow up in months
#'
#' @export
#'
stb_tl_simu_enroll <- function(ntot,
                               enroll_dur_mth,
                               min_fu      = NULL,
                               date_bos    = NULL,
                               mth_to_days = 30.4,
                               ...) {
    rand_enroll <- runif(ntot, 0, enroll_dur_mth)
    day_enroll  <- floor(rand_enroll * mth_to_days)

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

#' Simulate from piecewise constant exponential
#'
#'
#' @examples
#' rd_sim_pwexp(hazards = c("0.5" = 0.045, "Inf" = 0.025), offset = 0)
#'
#' @export
#'
stb_tl_surv_simu_pwexp <- function(hazards, offset = 0) {
    if (1 == length(hazards)) {
        tte <- rexp(1,  hazards)
    } else {
        segments <- as.numeric(names(hazards))
        segments <- sort(segments) - offset

        inx    <- min(which(segments > 0))
        cur_t  <- 0
        flag   <- FALSE
        while (!flag) {
            cur_h   <- hazards[inx]
            cur_int <- rexp(1, cur_h)

            if ((cur_t + cur_int) <= segments[inx]) {
                flag <- TRUE
                tte  <- cur_t + cur_int
                break
            }

            cur_int <- segments[inx]
            inx     <- inx + 1
        }
    }

    ## return
    rst <- c(tte    = tte,
             offset = offset)

    return(rst)
}

#' Log-Rank Test
#'
#' @return pvalue one-sided p-value
#'
#' @export
#'
stb_tl_surv_logrank <- function(data,
                                fml_survdiff =
                                    "Surv(day_pfs, status_pfs) ~ arm") {
    surv_diff <- coxph(as.formula(fml_survdiff), data = data)
    surv_sum  <- summary(surv_diff)
    pvalue    <- surv_sum$coefficients[1, 5]
    zscore    <- surv_sum$coefficients[1, 4]
    hr        <- surv_sum$coefficients[1, 2]

    c(zscore = zscore,
      pvalue = pnorm(zscore),
      nevent = surv_diff$nevent,
      hr     = hr)
}


#' Get Interim Analysis Data
#'
#' @param info_frac information fraction
#' @param total_events total number of events
#'
#' @export
#'
stb_tl_interim_data_2arm <- function(data,
                                     total_events,
                                     info_frac,
                                     event = "os") {
    v_date       <- paste("date_",   event, sep = "")
    v_status     <- paste("status_", event, sep = "")
    target_event <- floor(total_events * info_frac)

    ## all events
    events       <- data %>%
        dplyr::filter((!!sym(v_status)) == 1) %>%
        arrange(v_date)

    stopifnot(nrow(events) > target_event)

    ## date interim based on information fraction
    data$date_interim <- events[target_event, v_date]

    ## censor at interim
    rst <- data %>%
        filter(date_enroll  <= date_interim) %>%
        mutate(status_os = if_else(date_os <= date_interim,
                                   status_os,
                                   0),
               status_pfs = if_else(date_pfs <= date_interim,
                                    status_pfs,
                                    0),
               date_os = if_else(date_os <= date_interim,
                                 date_os,
                                 date_interim),
               date_pfs = if_else(date_pfs <= date_interim,
                                  date_pfs,
                                  date_interim),
               day_pfs = date_pfs - date_enroll,
               day_os  = date_os - date_enroll
               )

    rst
}
