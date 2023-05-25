## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for simulating recurrent events
##      All functions start with "stb_rc_"
##
##  DATE:
##      FEBRUARY, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#' @export
#'
rcurrent_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Study with recurrent events \n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(0.5, 0.5))\n")
    cat("    hr_by_arm:      hazard rate, i.e. annualized event rate, for each arm\n")
    cat("                    (default c(0.7, 0.3))\n")
    cat("    k_by_arm:       dispersion for each arm, smaller the bigger variance\n")
    cat("                    (default c(1, 1))\n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    max_n_event:    Maximum number of events (default 10)\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}


#' Default design parameter for single arm
#'
#'
internal_rcurrent_dpara <- function() {
    list(sample_size       = 100,
         ratio_by_arm      = c(0.5, 0.5),
         hr_by_arm         = c(1.3, 1),
         k_by_arm          = c(1, 1),
         par_enroll        = list(type = "by_duration", pt_dur_mth = 16),
         max_n_event       = 10,
         alpha             = 0.025
         )
}

#' Generate data
#'
#'
#'
rcurrent_gen_data <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n by arm
    n_by_arm <- tl_draw_arm_size(sample_size  = lst_design$sample_size,
                                 ratio_by_arm = lst_design$ratio_by_arm)

    ## enrollment
    dat_enroll <- do.call(stb_tl_simu_enroll_arms,
                          list(n_by_arm   = c(n_by_arm),
                               par_enroll = lst_design$par_enroll))

    ## outcome
    ## annual to day hr
    hr_by_arm   <- lst_design$hr_by_arm / 365.25
    k_by_arm    <- lst_design$k_by_arm
    max_n_event <- lst_design$max_n_event
    rst         <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_events <- stb_tl_rc_simu_events(n_by_arm[i],
                                            hr_by_arm[i],
                                            k_by_arm[i],
                                            n_event = max_n_event)

        cur_events$arm <- i - 1
        rst            <- rbind(rst, cur_events)
    }

    ## merget data
    rst <- dat_enroll %>%
        left_join(rst, by = c("arm", "sid" = "id")) %>%
        mutate(day_start = day_enroll + start,
               day_end   = day_enroll + end,
               sid       = sid + arm * lst_design$sample_size)


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
rcurrent_day_eos_1 <- function(data_full, min_fu_days = 12 * 7, pt_proportion = 1) {
    dat <- data_full %>%
        select(arm, sid, day_enroll) %>%
        distinct() %>%
        arrange(day_enroll)

    day_eos <- dat[ceiling(nrow(dat) * pt_proportion), "day_enroll"]
    day_eos <- day_eos + min_fu_days

    data_full %>%
        mutate(day_eos = day_eos)
}

#' Generate analysis data set
#'
#' Generate analysis dataset for fixed fu days
#'
rcurrent_day_eos_2 <- function(data_full, fu_days = 12 * 7, pt_proportion = 1) {
    dat <- data_full %>%
        select(arm, sid, day_enroll) %>%
        distinct() %>%
        mutate(day_eos = day_enroll + fu_days) %>%
        arrange(day_enroll)

    n_pt <- ceiling(nrow(dat) * pt_proportion)
    dat  <- dat[1 : n_pt, ]

    dat %>%
        left_join(data_full %>%
                  select(-day_enroll),
                  by = c("arm", "sid"))
}

#' Generate censoring data
#'
#'
rcurrent_censor <- function(data_full) {
    data_full %>%
        filter(day_start <= day_eos) %>%
        mutate(censor  = if_else(day_end <= day_eos, 0, 1),
               day_end = if_else(day_end <= day_eos, day_end, day_eos),
               time    = day_end - day_start)
}

#' Generate NB data
#'
#'
rcurrent_get_nb <- function(data_full) {
    n_arm <- max(data_full$arm)
    rst   <- NULL
    for (i in 0:n_arm) {
        cur_arm   <- data_full %>%
            filter(i == arm)

        cur_count <- cur_arm %>%
            group_by(arm, sid) %>%
            summarize(y = n() - 1)

        cur_rst <- cur_count %>%
            left_join(data_full %>%
                      select(sid, arm, lambda, day_enroll) %>%
                      distinct(),
                      by = c("arm", "sid"))

        rst <- rbind(rst, cur_rst)
    }

    rst
}
