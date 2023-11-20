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
    cat("    power:          power level (default 0.95)\n")
    cat("    max_n_event:    Maximum number of events for each patient (default 10)\n")
    cat("    par_enroll:     list of enrollment parameters \n")
    cat("    n_stage1:       no. of patients for stage 1 (default 100) \n")
    cat("    fix_fu:         fixed FU days (default 12 * 7) \n")
    cat("    rcur_info:      Information ratio in recurrent events (default 0.2) \n")
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
         alpha             = 0.025,
         power             = 0.9,
         n_stage1          = 100,
         fix_fu            = 12 * 7,
         rcur_info         = 0.2
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
rcurrent_day_eos_1 <- function(data_full,
                               min_fu_days   = 12 * 7,
                               pt_proportion = 1) {
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

#' Generate analysis data set
#'
#' Generate final analysis dataset for 1) stage 1 patients with fixed fu, 2)
#' stage 2 patients with fixed fu, 3) target number of events
#'
#'
rcurrent_day_eos_adapt_1 <- function(data_full, n_stage1, fix_fu = 12 * 7) {

    dat <- data_full %>%
        select(arm, sid, date_enroll) %>%
        distinct() %>%
        arrange(date_enroll)

    date_eos_1 <- dat[n_stage1, "date_enroll"] + 1
    rst        <- data_full %>%
        mutate(date_eos_1 = date_eos_1,
               date_eos_2 = date_enroll + fix_fu,
               date_eos   = if_else(date_eos_1 < date_eos_2,
                                    date_eos_1,
                                    date_eos_2),
               day_eos    = as.numeric(date_eos - date_bos)) %>%
        rcurrent_censor()

    rst
}


#' Generate analysis data set
#'
#' Generate final analysis dataset for 1) stage 1 patients with fixed fu, 2)
#' stage 2 patients with fixed fu, 3) target number of events
#'
#'
rcurrent_day_eos_adapt_2 <- function(data_full,
                                     n_stage1,
                                     n_stage2,
                                     target_event,
                                     rcur_info = 0.2,
                                     fix_fu    = 12 * 7) {

    dat <- data_full %>%
        select(arm, sid, date_enroll) %>%
        distinct() %>%
        arrange(date_enroll)

    dat_stage1 <- dat %>%
        slice(1:n_stage1) %>%
        left_join(data_full, by = c("arm" = "arm", "sid" = "sid"))

    dat_stage2 <- dat %>%
        slice(n_stage1 + 1:n_stage2) %>%
        left_join(data_full, by = c("arm" = "arm", "sid" = "sid"))

    ## last patient FU finished
    date_eos_1 <- dat_stage2[n_stage1, "date_enroll"] + fix_fu

    ## target event observed
    dat_stage2_target <- dat_stage2 %>%
        arrange(day_end) %>%
        mutate(nevent   = if_else(1 == inx, 1, rcur_info),
               cumevent = cumsum(nevent)) %>%
        filter(cumevent <= target_event) %>%
        slice_tail(1) %>%
        mutate(date_eos = date_bos + day_start + 1)

    date_eos_2 <- dat_stage2_target[1, "date_eos"]

    ## make sure stage 1 patients have FU
    date_eos_12 <- max(date_eos_1, date_eos_2)

    ##
    rst <- rbind(dat_stage1, dat_stage2) %>%
        mutate(date_eos_1 = date_eos_12,
               date_eos_2 = date_enroll + fix_fu,
               date_eos   = if_else(date_eos_1 < date_eos_2,
                                    date_eos_1,
                                    date_eos_2),
               day_eos    = as.numeric(date_eos - date_bos)) %>%
        rcurrent_censor()
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
                      select(sid, arm, lambda, day_enroll, day_eos) %>%
                      distinct(),
                      by = c("arm", "sid")) %>%
            mutate(day_onstudy = day_eos - day_enroll)

        rst <- rbind(rst, cur_rst)
    }

    rst
}
