## -----------------------------------------------------------------------------
##
## This file contains R functions for specifying a clinical study design.
##
## All functions start with "stb_design"
##
##
##
## -----------------------------------------------------------------------------

#' Get Study Design for Joint Primary and Secondary Endpoints
#'
#' @param btype_primary type for primary bound
#' @param bound_secondary alpha level at each interim analysis. Default is the
#'     naive alpha boundary.
#'
#'
#' @export
#'
stb_design_surv_join <- function(info_frac       = c(0.3, 0.7, 1),
                                 primary         = c("os", "pfs"),
                                 target_primary  = 100,
                                 sample_size     = 500,
                                 ctl_ratio       = 0.5,
                                 annual_drop     = 0.000001,
                                 enroll_dur_mth  = 18,
                                 alpha           = 0.025,
                                 btype_primary   = "asOF",
                                 bound_second    = NULL,
                                 ctl_median_os   = 15,
                                 ctl_median_pfs  = 9,
                                 hr_os           = 0.7,
                                 hr_pfs          = 1,
                                 rho_ctl         = 0.8,
                                 rho_trt         = NULL,
                                 method          = c("given.pfs", "given.os"),
                                 ...,
                                 seed            = NULL) {

    primary <- match.arg(primary)

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)


    ## hazard and correlation
    par_ctl <- stb_surv_join_par(median_os  = ctl_median_os,
                                 median_pfs = ctl_median_pfs,
                                 rho        = rho_ctl,
                                 method     = method,
                                 ...)

    if (is.null(rho_trt))
        rho_trt <- rho_ctl

    par_trt <- stb_surv_join_par(median_os  = ctl_median_os  / hr_os,
                                 median_pfs = ctl_median_pfs / hr_pfs,
                                 rho        = rho_trt,
                                 method     = method,
                                 ...)

    ## boundary
    bound_primary <-
        getDesignGroupSequential(sided            = 1,
                                 alpha            = alpha,
                                 informationRates = info_frac,
                                 typeOfDesign     = btype_primary)$stageLevels

    if (is.null(bound_second))
        bound_second  <- rep(alpha, length(bound_primary))

    pval_bounds           <- cbind(bound_primary, bound_second)
    colnames(pval_bounds) <- unique(c(primary, c("os", "pfs")))

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst         <- as.list(environment())
    class(rst)  <- "DESIGN_SURV_JOIN"

    rst
}


#' Get Study Design for Stratified Survival Analysis
#'
#' @param ctl_median median survival time for each stratum
#' @param hr hazard ratio for each stratum
#' @param stra_freq frequency of each stratum
#'
#' @export
#'
stb_design_surv_stra <- function(target          = 100,
                                 sample_size     = 500,
                                 ctl_ratio       = 0.5,
                                 annual_drop     = 0.000001,
                                 enroll_dur_mth  = 18,
                                 stra_freq       = c(0.5, 0.5),
                                 ctl_median      = c(15, 15),
                                 hr              = c(0.7, 0.7)) {

    ## return
    rst         <- as.list(environment())
    class(rst)  <- "DESIGN_SURV_STRATA"

    rst
 }


#' Get Study Design for PFS with Biomarker
#'
#' @param interim_biom Interim sample size for arm selection
#' @param f_arm_sel Arm selection function
#'
#' @export
#'
stb_design_surv_biom <- function(sample_size     = 600,
                                 ratio_by_arm    = c(1/3, 1/3, 1/3),
                                 p_biom_by_arm   = rbin(c(0.3, 0.7),
                                                        c(0.3, 0.7),
                                                        c(0.3, 0.7)),
                                 ctl_median_surv = c(10, 15),
                                 hr              = c(0.7, 0.7),
                                 annual_drop     = 0.000001,
                                 enroll_dur_mth  = 18,
                                 info_frac       = c(0.5, 0.7, 1),
                                 btype_primary   = "asOF",
                                 alpha           = 0.025,
                                 target_events   = 200,
                                 interim_biom    = 300,
                                 f_arm_sel       = stb_surv_biom_arm_sel_rule_1,
                                 fml_surv        = "Surv(day_pfs, status_pfs) ~ arm"
                                 ) {

    ## boundary
    bound <-
        getDesignGroupSequential(sided            = 1,
                                 alpha            = alpha,
                                 informationRates = info_frac,
                                 typeOfDesign     = btype_primary)$stageLevels

    ## median survival
    median_surv_by_arm <- ctl_median_surv
    for (n_arm in seq_len(length(hr))) {
        median_surv_by_arm <- rbind(
            median_surv_by_arm,
            ctl_median_surv / hr[n_arm])
    }

    ## return
    rst         <- as.list(environment())
    class(rst)  <- "DESIGN_SURV_BIOM"

    rst
}
