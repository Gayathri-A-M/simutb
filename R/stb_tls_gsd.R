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

    z_alpha <- qnorm(1 - alpha)
    theta   <- z_alpha + qnorm(power)
    b_t     <- sqrt(info_frac) * zscore

    if (use_observed)
        theta <- b_t / info_frac

    eb_1 <- b_t + theta * (1 - info_frac)
    rst  <- (z_alpha - eb_1) / sqrt(1 - info_frac)

    1 - pnorm(rst)
}


#' Compute predictive power
#'
#' Formula is from Proschan, Lan and Wittes (2006) Chapter 3
#'
#' @export
#'
stb_tl_gsd_predpower <- function(zscore, info_frac,
                                 alpha        = 0.025,
                                 power        = 0.9,
                                 prior_weight = 0.05) {

    stopifnot(prior_weight >= 0 &
              prior_weight <= 1)

    z_alpha <- qnorm(1 - alpha)
    theta   <- z_alpha + qnorm(power)
    b_t     <- sqrt(info_frac) * zscore

    sig0_2  <- (1 - prior_weight) / prior_weight

    rst     <- (b_t - z_alpha) * (1 + info_frac * sig0_2)
    rst     <- rst + (1 - info_frac) * (theta + b_t * sig0_2)
    denom   <- (1 - info_frac) * (1 + sig0_2) * (1 + info_frac * sig0_2)

    rst     <- rst / sqrt(denom)

    pnorm(rst)
}


#' GSD boundaries
#'
#'
#' @export
#'
stb_tl_gsd_boundary <- function(...,
                                info_fracs     = c(0.5),
                                alpha_spending = c(0.001),
                                sided          = 1,
                                alpha          = 0.025,
                                power          = 0.9) {

    info_fracs <- sort(info_fracs)
    if (max(info_fracs) < 1) {
        info_fracs     <- c(info_fracs, 1)
        alpha_spending <- c(alpha_spending, alpha)
    }

    stopifnot(length(info_fracs) == length(alpha_spending))

    rst <- getDesignGroupSequential(...,
                                    sided             = sided,
                                    alpha             = alpha,
                                    informationRates  = info_fracs,
                                    userAlphaSpending = alpha_spending)

    rst_rej <- stb_tl_gsd_prob(rst$informationRates,
                               boundary = rst$criticalValues,
                               boundary_type = "zscore",
                               alpha         = alpha,
                               power         = power)

    data.frame(inx                 = seq_len(length(rst$informationRates)),
               info_frac           = rst$informationRates,
               boundary_zscore     = rst$criticalValues,
               nominal_alpha       = rst$stageLevels,
               nominal_alpha_spent = alpha - rst$stageLevels,
               ia_alpha_spent      = diff(c(0, rst$alphaSpent)),
               cumu_alpha_spent    = rst$alphaSpent,
               ia_power            = rst_rej$ia_power,
               cumu_power          = rst_rej$cumu_power,
               study_alpha         = alpha,
               study_power         = power)
}


#' GSD covariance matrix
#'
#'
#' @export
#'
stb_tl_gsd_cov <- function(info_fracs = c(0.5, 1)) {

    n_ana   <- length(info_fracs)
    mat_sig <- matrix(NA, n_ana, n_ana)

    for (i in 1:n_ana) {
        for (j in i:n_ana) {
            mat_sig[i, j] <-
                mat_sig[j, i] <-
                sqrt(info_fracs[i] / info_fracs[j])
        }
    }
    diag(mat_sig) <- 1

    ## return
    mat_sig
}

#' GSD mean
#'
#'
#' @export
#'
stb_tl_gsd_mean <- function(info_fracs = c(0.5, 1),
                            alpha      = 0.025,
                            power      = 0) {
    if (power > 0) {
        z_alpha <- qnorm(1 - alpha)
        theta   <- z_alpha + qnorm(power)
    } else {
        theta <- 0
    }

    zscore_t <- theta * sqrt(info_fracs)

    ## return
    zscore_t
}

#' GSD boundaries
#'
#'
#' @export
#'
stb_tl_gsd_prob <- function(info_fracs    = c(0.5, 1),
                            boundary      = c(0.022, 0.022),
                            boundary_type = c("alpha",
                                              "alpha_spent",
                                              "zscore"),
                            alpha         = 0.025,
                            power         = 0.9) {

    boundary_type <- match.arg(boundary_type)
    stopifnot(max(info_fracs) == 1)
    stopifnot(length(info_fracs) == length(boundary))

    if ("zscore" == boundary_type) {
        boundary_zscore     <- boundary
        nominal_alpha       <- 1 - pnorm(boundary_zscore)
        nominal_alpha_spent <- alpha - nominal_alpha
    } else if ("alpha" == boundary_type) {
        nominal_alpha       <- boundary
        boundary_zscore     <- qnorm(1 - nominal_alpha)
        nominal_alpha_spent <- alpha - nominal_alpha
    } else {
        nominal_alpha_spent <- boundary
        nominal_alpha       <- alpha - nominal_alpha_spent
        boundary_zscore     <- qnorm(1 - nominal_alpha)
    }

    zscore_t <- stb_tl_gsd_mean(info_fracs,
                                alpha = alpha,
                                power = power)

    rst <- NULL
    for (i in seq_len(length(info_fracs))) {
        cur_mean <- zscore_t[1 : i]
        cur_sig  <- stb_tl_gsd_cov(info_fracs[1 : i])

        if (1 == i) {
            lb <- boundary_zscore[i]
            ub <- Inf
        } else {
            lb <- c(rep(-Inf, i - 1), boundary_zscore[i])
            ub <- c(boundary_zscore[1 : (i - 1)], Inf)
        }

        cur_rej_alpha <- pmvnorm(lower = lb,
                                 upper = ub,
                                 mean  = rep(0, i),
                                 sigma = cur_sig)

        if (0 == power) {
            cur_rej_pow <- cur_rej_alpha
        } else {
            cur_rej_pow <- pmvnorm(lower = lb,
                                   upper = ub,
                                   mean  = cur_mean,
                                   sigma = cur_sig)
}

        rst <- rbind(rst,
                     c(cur_rej_alpha,
                       cur_rej_pow))
    }


    ## return
    data.frame(
        inx                 = seq_len(length(info_fracs)),
        info_frac           = info_fracs,
        boundary_zscore     = boundary_zscore,
        nominal_alpha       = nominal_alpha,
        nominal_alpha_spent = nominal_alpha_spent,
        ia_alpha_spent      = rst[, 1],
        cumu_alpha_spent    = cumsum(rst[, 1]),
        ia_power            = rst[, 2],
        cumu_power          = cumsum(rst[, 2]),
        study_alpha         = alpha,
        study_power         = power)
}


#' GSD for one interim analysis
#'
#'
#' @export
#'
stb_tl_gsd_solve <- function(info_fracs = c(0.2, 1),
                             boundary   = c(0.001, NA),
                             boundary_type = c("alpha",
                                               "alpha_spent",
                                               "zscore"),
                             alpha  = 0.025,
                             power  = 0.9,
                             n_grid = 10000,
                             tol    = 1e-6) {


    stopifnot(max(info_fracs) == 1)
    stopifnot(length(info_fracs) == length(boundary))
    boundary_type <- match.arg(boundary_type)


    inx_na <- which(is.na(boundary))
    stopifnot(1 == length(inx_na))

    if ("zscore" == boundary_type) {
        vec <- seq(0, 3.5, length.out = n_grid)
    } else {
        vec <- seq(0, alpha, length.out = n_grid)
    }

    check <- parallel::mclapply(
        vec,
        function(x) {
            bd         <- boundary
            bd[inx_na] <- x
            cur_rst    <- stb_tl_gsd_prob(info_fracs    = info_fracs,
                                          boundary      = bd,
                                          boundary_type = boundary_type,
                                          alpha         = alpha,
                                          power         = 0)

            max(cur_rst$cumu_alpha_spent)
        },
        mc.cores = detectCores() - 1)

    check <- simplify2array(check)
    check <- abs(check - alpha)
    if (min(check) > tol) {
        return(NULL)
    }

    inx_min    <- which.min(check)
    bd         <- boundary
    bd[inx_na] <- vec[inx_min]
    rst        <- stb_tl_gsd_prob(info_fracs    = info_fracs,
                                  boundary      = bd,
                                  boundary_type = boundary_type,
                                  alpha         = alpha,
                                  power         = power)


    rst$inx_na <- inx_na

    ## return
    rst
}
