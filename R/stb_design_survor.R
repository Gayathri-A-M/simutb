## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating oncology trials with OR status
##
##
##
## -----------------------------------------------------------------------------


#' Oncology Study Design I
#'
#'
survor_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Oncology single arm trial with OR status and survival\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    p_or:           Probability of OR (default 0.4)\n")
    cat("    or_median_prog: Median months to progression for patients\n")
    cat("                    with OR (default 6) \n")
    cat("    or_median_dth:  Median months to death for patients\n")
    cat("                    with OR (default 10) \n")
    cat("    nor_median_prog:Median months to progression for patients\n")
    cat("                    without OR (default 4) \n")
    cat("    nor_median_dth: Median months to death for patients\n")
    cat("                    without OR (default 6) \n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    min_fu:         minimum follow-up months (default 6) \n")
    cat("    enroll_dur_mth: enrollment months (default 18) \n")
    cat("    date_bos:       begin of study date\n")
    cat("                    (default 2022-1-1) \n")
}

#' Default design paramters
#'
#'
survor_default_para <- function() {
    lst_default <- list(sample_size      = 100,
                        enroll_dur_mth   = 18,
                        date_bos         = as.Date("2022-1-1"),
                        min_fu           = 6,
                        p_or             = 0.4,
                        or_median_prog   = 6,
                        or_median_dth    = 10,
                        nor_median_prog  = 4,
                        nor_median_dth   = 6,
                        annual_drop      = 0.00001)
}

#' Simulate a trial
#'
#' @param sample_size Total number of patients
#' @param enroll_dur_mth Duration of enrollment in months
#' @param min_fu Minimum follow up months
#' @param date_bos Begin of study date
#' @param mth_to_days Days in each month
#' @param p_or Probability of overall response
#' @param or_median_prog median month for progression for OR patients
#' @param or_median_dth median month for death for OR patients
#' @param nor_median_prog median month for progression for patients without OR
#' @param nor_median_dth median month for death for OR patients without OR
#' @param annual_drop annual drop rate
#'
#'
#' @export
#'
survor_gen_data <- function(sample_size      = 100,
                           enroll_dur_mth   = 18,
                           date_bos         = as.Date("2022-1-1"),
                           min_fu           = 6,
                           p_or             = 0.4,
                           or_median_prog   = 6,
                           or_median_dth    = 10,
                           nor_median_prog  = 4,
                           nor_median_dth   = 6,
                           annual_drop      = 0.00001,
                           ...) {

    ## enrollment
    dta_enroll <- stb_tl_simu_enroll_arc(sample_size, enroll_dur_mth,
                                         date_bos = date_bos, min_fu = min_fu)

    ## overall response
    dta_or     <- stb_simu_or(dta_enroll, p_or = p_or)

    ## patients with OR
    sub_or     <- stb_simu_event(dta_or %>% filter(status_or == TRUE),
                                 or_median_prog,
                                 or_median_dth,
                                 annual_drop,
                                 shift = "day_or")

    ## patients without OR
    sub_nor    <- stb_simu_event(dta_or %>% filter(status_or == FALSE),
                                 nor_median_prog,
                                 nor_median_dth,
                                 annual_drop,
                                 shift = NULL)

    rst <- rbind(sub_or, sub_nor)

    ## get date
    if (!is.null(date_bos)) {
        rst <- rst %>%
            mutate(date_prog   = date_enroll + day_prog,
                   date_dth    = date_enroll + day_dth,
                   date_censor = date_enroll + day_censor,
                   date_pfs    = date_enroll + day_pfs,
                   date_os     = date_enroll + day_os)
    }

    ## return
    rst
}


#' Simulate Overall Response
#'
#' @inheritParams survor_gen_data
#'
#' @export
#'
stb_simu_or <- function(dta_enroll, p_or) {
    rand_or <- rbinom(nrow(dta_enroll), 1, p_or)
    day_or  <- sapply(rand_or, function(i) {
        if (0 == i) {
            rst <- NA
        } else {
            rst <- runif(1,
                         min = 1,
                         max = dta_enroll$day_eos[i])
            rst <- floor(rst)
        }
    })

    dta_enroll$status_or <- 1 == rand_or
    dta_enroll$day_or    <- day_or + dta_enroll$day_enroll

    if (!is.null(dta_enroll$date_enroll)) {
        dta_enroll$date_or <- dta_enroll$day_or  + dta_enroll$date_enroll
    }

    dta_enroll
}


#' Simulate progression, death, and censoring
#'
#' @param ntot  total number of patients
#' @param median_prog median month for progression
#' @param median_dth median month for death
#' @param annual_drop annual drop rate
#'
#' @export
#'
stb_simu_event <- function(dta_or,
                           median_prog,
                           median_dth,
                           annual_drop,
                           shift = "day_or") {

    ntot <- nrow(dta_or)

    if (is.null(shift)) {
        day_shift <- 0
    } else {
        day_shift <- dta_or[[shift]]
    }

    ## time to prog
    day_prog   <- day_shift + stb_tl_rexp(ntot, median_mth  = median_prog)
    day_dth    <- day_shift + stb_tl_rexp(ntot, median_mth  = median_dth)
    day_censor <- day_shift + stb_tl_rexp(ntot, annual_drop = annual_drop)
    day_censor <- apply(cbind(dta_or$day_eos, day_censor), 1, min)

    ## pfs and os
    dta_pfs_os <- stb_tl_pfs_os(day_prog, day_dth, day_censor)

    cbind(dta_or, dta_pfs_os)
}

#' Simulate trials
#'
#' @param ... Parameter for f_single
#'
#' @export
#'
stb_simu_trials <- function(nreps = 100,
                            f_single = stb_simu_trial_single, ...,
                            n_cores = 5) {
    rst <- parallel::mclapply(1:nreps,
                              function(k) {
                                  f_single(...)
                              }, mc.cores = n_cores)

    rst
}
