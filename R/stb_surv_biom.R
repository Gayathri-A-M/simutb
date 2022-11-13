## -----------------------------------------------------------------------------
##
## This file contains R functions for simulating joint PFS and a biomarker.
##
## All functions start with "stb_surv_biom"
##
##
##
## -----------------------------------------------------------------------------


#' Simulate an arm
#'
#' @param p_biomarker biomarker rates
#' @param median_mth median survival for responders and non-responsders
#' @export
#'
stb_surv_biom_arm_simu <- function(n,
                                   p_biomarker,
                                   median_mth,
                                   enroll_dur_mth,
                                   annual_drop = 0,
                                   date_bos    = as.Date("2022-1-1"),
                                   seed        = NULL,
                                   ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    smp_biom   <- rmultinom(1, n, p_biomarker)
    dta_enroll <- stb_tl_simu_enroll(n,
                                     enroll_dur_mth,
                                     date_bos    = date_bos,
                                     ...)
    dta_censor <- stb_tl_rexp(n,
                              median_mth  = NULL,
                              hazard      = NULL,
                              annual_drop = annual_drop,
                              ...)

    dta_biom  <- NULL
    dta_event <- NULL
    for (j in seq_len(length(smp_biom))) {
        cur_biom  <- rep(j - 1, smp_biom[j])
        cur_event <- stb_tl_rexp(smp_biom[j],
                                 median_mth  = median_mth[j],
                                 hazard      = NULL,
                                 annual_drop = NULL,
                                 ...)

        dta_biom  <- c(dta_biom,  cur_biom)
        dta_event <- c(dta_event, cur_event)
    }

    days <- stb_tl_pfs_os(dta_event,
                          dta_event,
                          dta_censor)

    rst <- cbind(dta_enroll, days, biom = dta_biom) %>%
        data.frame() %>%
        select(-day_prog, -day_dth, -day_event) %>%
        mutate(date_censor = date_enroll + day_censor,
               date_pfs    = date_enroll + day_pfs,
               date_os     = date_enroll + day_os)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Simulate a trial
#'
#'
#' @export
#'
stb_surv_biom_trial_simu <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n_by_arm    <- rmultinom(1,
    ##                          lst_design$sample_size,
    ##                          lst_design$ratio_by_arm)
    n_by_arm           <- lst_design$sample_size
    n_by_arm           <- n_by_arm * lst_design$ratio_by_arm
    n_by_arm           <- n_by_arm / sum(lst_design$ratio_by_arm)
    n_by_arm           <- floor(n_by_arm)

    p_biom_by_arm      <- lst_design$p_biom_by_arm
    median_surv_by_arm <- lst_design$median_surv_by_arm

    rst <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_arm <- stb_surv_biom_arm_simu(
            n              = n_by_arm[i],
            p_biomarker    = p_biom_by_arm[i, ],
            median_mth     = median_surv_by_arm[i, ],
            enroll_dur_mth = lst_design$enroll_dur_mth,
            annual_drop    = lst_design$annual_drop,
            ...)

        cur_arm$arm <- as.character(i - 1)
        rst         <- rbind(rst, cur_arm)
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Plot a simulated data
#'
#'
#' @export
#'
stb_surv_biom_trial_plot <- function(data, facet_by = "arm") {
    fit      <- survfit(Surv(day_pfs, status_pfs) ~ arm + biom, data = data)
    plt_pfs  <- ggsurvplot_facet(fit,
                                 data = data,
                                 facet.by = c(facet_by))

    plt_pfs
}

#' Arm Selection Rule I
#'
#' Arm selection by interim biomarker. Select the treatment arm with the largest
#' biomarker mean.
#'
#' @export
#'
stb_surv_biom_arm_sel_rule_1 <- function(data) {
    dta_biom <- data %>%
        filter(arm > 0 & included == 1) %>%
        group_by(arm) %>%
        summarize(n         = n(),
                  mean_biom = mean(biom)) %>%
        arrange(mean_biom)

    sel_arm <- as.numeric(dta_biom[nrow(dta_biom),
                                   'arm'])

    data    <- data %>% filter(arm %in% as.character(c(0, sel_arm)))

    list(sel_arm = sel_arm,
         biom    = dta_biom,
         data    = data)
}

#' Simulate a trial with interim
#'
#'
#'
#' @export
#'
stb_surv_biom_trial_interim <- function(lst_design, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    data_full <- stb_surv_biom_trial_simu(lst_design)

    ## check
    ## data_full %>%
    ##     group_by(arm, biom) %>%
    ##     summarize(n = n(),
    ##               pfs = median(day_pfs) / 30.1)

    ## arm selection
    f_sel        <- lst_design$f_arm_sel
    interim_date <- data_full %>% arrange(date_enroll)
    interim_date <- interim_date[lst_design$interim_biom,
                                 "date_enroll"]

    data_biom    <- data_full %>%
        mutate(included = if_else(date_enroll <= interim_date,
                                  1, 0))
    lst_arm_sel <- f_sel(data_biom)

    ## interim analysis
    info_frac <- lst_design$info_frac
    rst       <- NULL
    for (i in seq_len(length(lst_design$info_frac))) {
        data_interim <- stb_tl_interim_data_2arm(lst_arm_sel$data,
                                                 lst_design$target_events,
                                                 info_frac[i],
                                                 event = "pfs")

        test_pfs <- stb_tl_surv_logrank(data = data_interim,
                                        fml  = lst_design$fml_surv)

        pval     <- unname(test_pfs["pvalue_oneside"])
        bound    <- lst_design$bound[i]
        rst      <- rbind(rst,
                          c(arm        = lst_arm_sel$sel_arm,
                            inx        = i,
                            info_frac  = info_frac[i],
                            bound      = bound,
                            pval       = pval,
                            rej        = pval <= bound,
                            zscore     = unname(test_pfs["zscore"]),
                            nevent     = unname(test_pfs["nevent"]),
                            hr         = unname(test_pfs["hr"]))
                          )
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rownames(rst) <- NULL
    list(rst     = data.frame(rst),
         biom    = lst_arm_sel$biom)
}


#' Summarize simulation results
#'
#' @export
#'
stb_surv_biom_summary <- function(lst_rst) {

    n_reps      <- length(lst_rst)
    rst_biom    <- list()
    rst_interim <- list()
    for (i in 1:length(lst_rst)) {
        rst_biom[[i]]    <- lst_rst[[i]]$biom
        rst_interim[[i]] <- lst_rst[[i]]$rst %>% mutate(rep = i)
    }
    rst_biom    <- rbindlist(rst_biom)
    rst_interim <- rbindlist(rst_interim)

    ## arm selection
    rst_sel_arm <- rst_interim %>%
        distinct(arm, rep) %>%
        group_by(arm) %>%
        summarize(n    = n(),
                  rate = n / n_reps)

    ## biomarker at interim
    rst_biom_summary <- rst_biom %>%
        group_by(arm) %>%
        summarize(n         = mean(n),
                  mean_biom = mean(mean_biom))

    ## estimate
    rst_hr <- rst_interim %>%
        group_by(arm) %>%
        summarize(hr = mean(hr))

    ## rejection
    rej <- rst_interim %>%
        filter(rej == 1) %>%
        group_by(rep, arm) %>%
        arrange(inx) %>%
        slice(n = 1) %>%
        ungroup() %>%
        group_by(inx, info_frac, arm) %>%
        summarize(N_Rej = n()) %>%
        mutate(Rej = N_Rej / n_reps) %>%
        group_by(arm) %>%
        mutate(CumuRej = cumsum(Rej))

    rej_any <- rej %>%
        group_by(inx, info_frac) %>%
        summarize(N_Rej   = sum(N_Rej),
                  Rej     = sum(Rej),
                  CumuRej = sum(CumuRej))

    ## return
    rst <- list(sel_arm     = rst_sel_arm,
                biom        = rst_biom_summary,
                rej_arm     = rej,
                rej_any     = rej_any,
                hr          = rst_hr,
                raw_interim = rst_interim,
                raw_biom    = rst_biom)

    rst
}
