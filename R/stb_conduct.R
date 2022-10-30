## -----------------------------------------------------------------------------
##
## This file contains R functions for specifying conducting a clinical study
## design in simulation studies.
##
## All functions start with "stb_conduct"
##
##
##
## -----------------------------------------------------------------------------

#' Conduct Simulation Study for Joint Primary and Secondary Endpoints
#'
#'
#' @export
#'
stb_conduct_surv_join <- function(lst_design,
                                  n_rep  = 1000,
                                  n_core = 1,
                                  seed = NULL) {

    f_each <- function(k) {
        rst_interim <- stb_surv_join_trial_interim(
            ntrt           = ceiling(lst_design$sample_size *
                                     (1 - lst_design$ctl_ratio)),
            par_trt        = lst_design$par_trt,
            nctl           = floor(lst_design$sample_size *
                                   lst_design$ctl_ratio),
            par_ctl        = lst_design$par_ctl,
            enroll_dur_mth = lst_design$enroll_dur_mth,
            annual_drop    = lst_design$annual_drop,
            total_events   = lst_design$target_primary,
            info_frac      = lst_design$info_frac,
            primary        = lst_design$primary,
            seed           = all_seeds[k])

        rst_rej <- stb_surv_join_trial_hiertest(
            rst_interim,
            pval_bounds = lst_design$pval_bounds)

        rst_rej$rep <- k
        rst_rej
    }

    f_rej <- function(endp) {
        tmp <- rst_summary$rejection %>%
            filter(event == endp)

        if (0 == nrow(tmp)) {
            rst <- 0
        } else {
            rst <- max(tmp$CumuRej)
        }

        rst
    }

    stopifnot("DESIGN_SURV_JOIN" %in% class(lst_design))

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## all random seeds
    all_seeds <- ceiling(abs(rnorm(n_rep) * 10000))

    ## replications
    rst <- parallel::mclapply(seq_len(n_rep),
                              function(k) {
                                  if (0 == k %% 50)
                                      print(k)

                                  f_each(k)
                              },
                              mc.cores = n_core)

    rst <- rbindlist(rst)

    ## summary
    rst_summary <- stb_surv_join_summary(rst,
                                         colnames(lst_design$pval_bounds))
    rst_key     <- data.frame(
        pri_sec        = paste(colnames(lst_design$pval_bounds),
                               collapse = ","),
        info_frac      = paste(lst_design$info_frac,
                               collapse = ","),
        hr_os          = lst_design$hr_os,
        hr_pfs         = lst_design$hr_pfs,
        sample_size    = lst_design$sample_size,
        target_primary = lst_design$target_primary,
        rho_ctl        = lst_design$rho_ctl,
        rho_trt        = lst_design$rho_trt,
        kendall_ctl    = lst_design$par_ctl$kendall,
        kendall_trt    = lst_design$par_trt$kendall,
        hz_prog_ctl    = lst_design$par_ctl$hazard_prog,
        hz_prog_trt    = lst_design$par_trt$hazard_prog,
        zscore_cor     = rst_summary$zscore_cor,
        rej_os         = f_rej("os"),
        rej_pfs        = f_rej("pfs"),
        seed           = seed)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    list(rst_raw     = rst,
         rst_summary = rst_summary,
         rst_key     = rst_key)
}
