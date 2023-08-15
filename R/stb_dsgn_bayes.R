## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      Bayesian Design
##
##  DATE:
##      FEBRUARY, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
bayes_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Bayesian single-arm or two-arm design for binary outcomes.\n\n")

    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 50)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm\n")
    cat("                    (default c(0.5, 0.5))\n")
    cat("    p_by_arm:       response rate for each arm\n")
    cat("                    (default c(0.2, 0.2))\n")
    cat("    prior_by_arm:   prior a and b by arm\n")
    cat("                    (default rbind(c(1, 1), c(1, 1))\n")
    cat("    n_post:         number of posterior samples \n")
    cat("                    (default 3000)\n")
    cat("    x_post:         vector of x for calculating posterior pdf \n")
    cat("                    (default seq(0, 1, by = 0.01))\n")
    cat("    decision_h0:    Threshold under null hypothesis \n")
    cat("                    (default 0.2)\n")
    cat("    decision_gl:    Threshold greater.than or less.than \n")
    cat("                    (default 'greater.than')\n")
    cat("    decision_thresh:Threshold of success \n")
    cat("                    (default 0.6)\n")
    cat("    decision_ref:   Reference arm \n")
    cat("                    (default 1)\n")
}


#' Default design parameter for single arm
#'
#'
internal_bayes1arm_dpara <- function() {
    list(sample_size     = 50,
         ratio_by_arm    = 1,
         p_by_arm        = 0.2,
         prior_by_arm    = rbind(c(1, 1)),
         n_post          = 3000,
         x_post          = seq(0, 1, by = 0.01),
         decision_h0     = 0.2,
         decision_gl     = "greater.than",
         decision_thresh = 0.8,
         decision_ref    = 1)
}

internal_bayes2arm_dpara <- function() {
    list(sample_size     = 50,
         ratio_by_arm    = c(0.5, 0.5),
         p_by_arm        = c(0.2, 0.2),
         prior_by_arm    = rbind(c(1, 1),
                                 c(1, 1)),
         n_post          = 3000,
         x_post          = seq(0, 1, by = 0.01),
         decision_h0     = 0.2,
         decision_gl     = "greater.than",
         decision_thresh = 0.8,
         decision_ref    = 1)
}

#' Generate data
#'
#'
#'
bayes_gen_data <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n by arm
    n_by_arm <- tl_draw_arm_size(sample_size  = lst_design$sample_size,
                                 ratio_by_arm = lst_design$ratio_by_arm)

    ## response
    rst <- NULL
    for (i in seq_len(length(n_by_arm))) {
        y   <- rbinom(1, n_by_arm[i], lst_design$p_by_arm[i])
        rst <- rbind(rst,
                     c(arm = i,
                       n   = n_by_arm[i],
                       y   = y,
                       p   = lst_design$p_by_arm[i]))
    }

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    data.frame(rst)
}

#' Posterior Distribution
#'
#'
#' @export
#'
bayes_ana_data <- function(data, prior_by_arm,  ...) {

    n_arm    <- nrow(data)
    rst_post <- list()
    for (i in seq_len(n_arm)) {
        rst_post[[i]] <- stb_tl_bayes_beta_post(obs_y  = data[i, "y"],
                                                obs_n  = data[i, "n"],
                                                pri_ab = prior_by_arm[i, ],
                                                ...)
    }

    rst_post
}

#' Decision making
#'
#'
#' @export
#'
bayes_ana_decision <- function(rst_post, decision_ref = 1, ...) {

    n_arm <- length(rst_post)

    ## difference
    rst_diff <- list()
    if (1 == n_arm) {
        rst_diff    <- list(rst_post[[1]]$post_smp)
    } else {
        post_ref <- rst_post[[decision_ref]]$post_smp
        for (i in seq_len(n_arm)) {
            rst_diff[[i]] <- rst_post[[i]]$post_smp - post_ref
        }
    }

    ## success
    rst_success <- NULL
    for (i in seq_len(n_arm)) {
        cur_rst       <- stb_tl_bayes_decision(rst_diff[[i]],
                                               ...)
        rst_success   <- rbind(rst_success, c(arm = i, cur_rst))
    }

    list(rst_diff    = rst_diff,
         rst_success = rst_success)
}
