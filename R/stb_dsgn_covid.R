## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##
##      This file contains R functions for multi-stage multi-arm study
##      with survival endpoint
##
##  DATE:
##
##      MAY, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#' @export
#'
covid_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Multi-Stage Multi-Arm Trial with Survival Outcome\n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(1/3, 1/3, 1/3))\n")
    cat("    six_mth_cumu:   six month cumulative incidence rate for\n")
    cat("                    each arm (default c(0.021, 0.005, 0.005))\n")
    cat("    annual_drop:    annual dropout rate (default 0.000001) \n")
    cat("    min_fix_fu:     fixed follow-up months (default 6) \n")
    cat("    target_primary: target primary event (default 54) \n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    power:          pairwise comparison power (default 0.9)\n")
    cat("    par_interim:    list of parameters for creating interim \n")
    cat("                    and final analysis dataset\n")
    cat("    par_analysis:   list of analysis parameters\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}
