## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for simulating enrollment time
##
##  DATE:
##      APRIL, 2023
## -----------------------------------------------------------------------------


#' Simulate Repeated Measure
#'
#' @export
#'
stb_tl_simu_rm <- function(ntot, mat_mu, mat_sig) {
    rmvnorm(ntot,
            mean  = mat_mu,
            sigma = mat_sig)
}
