## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for group sequential designs
##
##  DATE:
##      AUG 2023
## -----------------------------------------------------------------------------

#' Compute conditional power
#'
#' Formula is from Proschan, Lan and Wittes (2006) Chapter 3
#'
#' @export
#'
stb_tl_gsd_condpower <- function(zscore, info_frac,
                                 alpha        = 0.025,
                                 power        = 0.9,
                                 use_observed = TRUE) {

    z_alpha   <- qnorm(1 - alpha)
    theta     <- z_alpha + qnorm(power)
    b_t       <- sqrt(info_frac) * zscore

    if (use_observed)
        theta <- b_t / info_frac

    eb_1 <- b_t + theta * (1 - info_frac)
    rst  <- (z_alpha - eb_1) / sqrt(1 - info_frac)

    1 - pnorm(rst)
}
