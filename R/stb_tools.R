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
                        hazard       = NULL,
                        median_mth   = 5,
                        annual_drop  = NULL,
                        mth_to_days  = 30.4,
                        take_floor   = TRUE) {

    if (is.null(hazard)) {
        if (!is.null(median_mth)) {
            hazard <- - log(0.5)  / median_mth
        } else {
            hazard <- - log(1 - annual_drop) / 12
        }
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
stb_tl_pfs_os <- function(day_prog, day_dth, day_censor) {
    f_s <- function(prog, dth, censor) {
        day_pfs <- min(prog, dth, censor)
        day_os  <- min(dth, censor)

        status_pfs <- censor > day_pfs
        status_os  <- censor > day_os

        c(prog,
          dth,
          censor,
          day_pfs,
          status_pfs,
          day_os,
          status_os)
    }

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
                       "status_os")

    rst
}
