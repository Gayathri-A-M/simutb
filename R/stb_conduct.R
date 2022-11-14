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
            ntrt           = ntrt,
            par_trt        = lst_design$par_trt,
            nctl           = nctl,
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

    stopifnot("DESIGN_SURV_JOIN" %in% class(lst_design))

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## all random seeds
    all_seeds <- ceiling(abs(rnorm(n_rep) * 100000))

    ## arm size
    nctl <- floor(lst_design$sample_size * lst_design$ctl_ratio)
    ntrt <- lst_design$sample_size - nctl

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
    rst_summary <- stb_surv_join_summary(rst, colnames(lst_design$pval_bounds))
    ##_key
    rst_key     <- stb_surv_join_key(lst_design, rst_summary, seed)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    list(rst_raw     = rst,
         rst_summary = rst_summary,
         rst_key     = rst_key)
}



#' Conduct Simulation Study for Arm Selection
#'
#'
#' @export
#'
stb_conduct_surv_biom <- function(lst_design,
                                  n_rep  = 1000,
                                  n_core = 1,
                                  ...,
                                  seed   = NULL) {

    stopifnot("DESIGN_SURV_BIOM" %in% class(lst_design))

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## all random seeds
    all_seeds <- ceiling(abs(rnorm(n_rep) * 10000))

    ## replication
    lst_rst <- parallel::mclapply(seq_len(n_rep),
                                  function(i) {
                                      if (i %% 50 == 0)
                                          print(i)

                                      stb_surv_biom_trial_interim(
                                          lst_design,
                                          seed = all_seeds[i])
                                  },
                                  mc.cores = n_core)

    ## summary
    rst_summary <- stb_surv_biom_summary(lst_rst)

    ## key
    rst_key     <- stb_surv_biom_key(lst_design, rst_summary, seed)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    list(rst_summary = rst_summary,
         rst_key     = rst_key)
}
