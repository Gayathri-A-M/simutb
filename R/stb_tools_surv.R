## -----------------------------------------------------------------------------
##
## This file contains R tool functions for survival analysis.
##
## All functions start with "stb_tl_".
##
##
##
## -----------------------------------------------------------------------------

#' Bivariate normal
#'
#' @export
#'
stb_tl_binorm <- function(n, rho) {
    ## random bivariate normal samples
    rnd_smp_1 <- rnorm(n)
    rnd_smp_2 <- rho * rnd_smp_1 + sqrt(1 - rho^2) * rnorm(n)

    ## random cdf
    rnd_cdf_1  <- pnorm(rnd_smp_1)
    rnd_cdf_2  <- pnorm(rnd_smp_2)

    list(rnd_smp = cbind(rnd_smp_1, rnd_smp_2),
         rnd_cdf = cbind(rnd_cdf_1, rnd_cdf_2))
}

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
                        take_floor   = FALSE) {

    if (is.null(hazard)) {
        hazard <- stb_tl_hazard(median_surv = median_mth,
                                annual_drop = annual_drop)[1]
    }

    if (0 == hazard) {
        rand_event <- rep(Inf, ntot)
    } else {
        rand_event <- rexp(ntot, hazard) * mth_to_days
    }

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
        day_event  <- min(prog, dth)

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

    rst           <- t(rst)
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
    day_enroll  <- rand_enroll * mth_to_days

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
                                    "Surv(day_pfs, status_pfs) ~ arm",
                                method = c("logrank", "score")) {

    method    <- match.arg(method)
    fml       <- as.formula(fml_survdiff)

    surv_diff <- survdiff(fml, data = data)
    pval_lr   <- surv_diff$pvalue
    surv_diff <- coxph(fml, data = data)
    surv_sum  <- summary(surv_diff)
    pval_cox  <- unname(surv_sum$sctest[3])
    hr        <- surv_sum$coefficients[, 2]
    zscore    <- surv_sum$coefficients[, 4]
    pval_z    <- pnorm(zscore)

    inx           <- seq_len(length(hr))
    inx[1]        <- ""
    names(hr)     <- paste("hr",     inx, sep = "")
    names(zscore) <- paste("zscore", inx, sep = "")
    names(pval_z) <- paste("pval_z_ones", inx, sep = "")

    ## one-sided pvalue
    pvalue <- switch(method,
                     logrank = pval_lr,
                     score   = pval_cox)

    if (hr[1] < 1) {
        pvalue <- pvalue / 2
    } else {
        pvalue <- 1 - pvalue / 2
    }

    ## return
    c(hr,
      zscore,
      pval_z,
      pvalue_oneside = pvalue,
      nevent         = surv_diff$nevent,
      pval_lr        = pval_lr,
      pval_cox       = pval_cox)
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
        arrange(!!sym(v_date))

    stopifnot(nrow(events) >= target_event)

    ## date interim based on information fraction
    data$date_interim <- events[target_event, v_date]

    ## censor at interim
    rst <- data %>%
        filter(date_enroll <= date_interim) %>%
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
               day_os  = date_os  - date_enroll
               )
    rst
}

#' Get Interim Analysis Data
#'
#' This is an improved version of stb_tl_interim_data that allows to create an
#' interim snapshot based on number of events or samples from any specific arms
#'
#' @param info_frac information fraction
#' @param total total number of events or samples
#' @param offset_days offset days after the target event for creating the
#'     interim data
#'
#' @export
#'
stb_tl_interim_data <- function(data, info_frac,
                                total       = NULL,
                                event       = c("os", "pfs", "enroll"),
                                arms        = NULL,
                                offset_days = 0) {

    event  <- match.arg(event)
    v_date <- paste("date_",   event, sep = "")

    ## default total is the total sample size
    if (is.null(total))
        total <- nrow(data)

    ## target number of events or samples
    target <- floor(total * info_frac)

    ## all events
    dta_events <- data
    if (event %in% c("os", "pfs")) {
        v_status   <- paste("status_", event, sep = "")
        dta_events <- dta_events %>%
            dplyr::filter((!!sym(v_status)) == 1)
    }

    if (!is.null(arms)) {
        dta_events <- dat_events %>%
            dplyr::filter(arm %in% arms)
    }

    dta_events  <- dta_events %>%
        arrange(!!sym(v_date))

    stopifnot(nrow(dta_events) >= target)

    ## censor at interim
    rst <- data %>%
        mutate(date_interim =
                   dta_events[target, v_date] +
                   offset_days) %>%
        filter(date_enroll <= date_interim) %>%
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
               day_os  = date_os  - date_enroll
               ) %>%
        arrange(date_enroll)
    rst
}
