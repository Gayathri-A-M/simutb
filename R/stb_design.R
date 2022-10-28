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
                                 ...,
                                 seed            = NULL) {

    primary <- match.arg(primary)

    ## seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)


    ## hazard and correlation
    par_ctl <- stb_surv_join_par(median_os  = ctl_median_os,
                                 median_pfs = ctl_median_pfs,
                                 rho        = rho_ctl, ...)

    if (is.null(rho_trt))
        rho_trt <- rho_ctl

    par_trt <- stb_surv_join_par(median_os  = ctl_median_os  / hr_os,
                                 median_pfs = ctl_median_pfs / hr_pfs,
                                 rho        = rho_trt,
                                 ...)

    ## boundary
    bound_primary <-
        getDesignGroupSequential(sided            = 1,
                                 alpha            = alpha,
                                 informationRates = info_frac,
                                 typeOfDesign     = btype_primary)$stageLevels

    if (is.null(bound_second))
        bound_second  <- rep(alpha, length(bound_primary))

    pval_bounds <- cbind(bound_primary, bound_second)

    colnames(pval_bounds) <- unique(c(primary, c("os", "pfs")))

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst         <- as.list(environment())
    class(rst)  <- "DESIGN_SURV_JOIN"

    rst
}
