## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for simulating recurrent events
##      All functions start with "stb_rc_"
##
##  DATE:
##      FEBRUARY, 2023
## -----------------------------------------------------------------------------

#' Power calculation
#'
#' Power calculation using the formula proposed in Zhu and Lakkis (2014)
#'
#' @param n total sample size, assuming 1:1 randomization ratio
#' @param mu_t mean follow-up time
#' @param r0 mean hazard rate in control arm
#' @param r1 mean hazard rate in treatment arm
#' @param k dispersion parameter in negative binomial
#' @param alpha type I error rate
#'
#' @examples
#'    stb_rc_power(n = 100, mu_t = 2, r0 = 3, r1 = 1.9, k = 1)
#'
#' @export
#'
stb_rc_power <- function(n, mu_t, r0, r1, k = 1, alpha = 0.05) {
    v_0 <- 2 / mu_t / r0 + 2 * k
    v_1 <- (1 / r0 + 1 / r1) / mu_t + 2 * k

    ss  <- sqrt(n / 2) * log(r0 / r1)
    ss  <- ss - sqrt(v_0) * qnorm(1 - alpha / 2)
    ss  <- ss / sqrt(v_1)

    pnorm(ss)
}


#' Minimal detectable difference
#'
#' Inverse power calculation to find the minimal detectable r1 given the rest of
#' the parameters
#'
#' @param target_power target power level. Default 0.5 for minimal detectable
#'     difference
#' @param var_solve parameter to be solved
#' @param ... Parameters for \link{stb_rc_power}
#'
#'
#' @examples
#'    r1 <- stb_rc_mdd(n = 100, mu_t = 2, r0 = 3, k = 1)
#'
#' @export
#'
stb_rc_mdd <- function(var_solve    = c("r1", "r0", "n", "mu_t"),
                       target_power = 0.5,
                       ...) {

    ft <- function(x) {
        pars_x[[var_solve]] <- x
        cur_pars            <- c(pars, pars_x)
        cur_power           <- do.call(stb_rc_power, cur_pars)

        cur_power - target_power
    }

    var_solve <- match.arg(var_solve)
    pars      <- list(...)
    pars_x    <- list()

    ## range
    var_range <- switch(var_solve,
                        r1   = c(1e-6, pars$r0),
                        r0   = c(pars$r1, 1e6),
                        n    = c(1e-6, 1e6),
                        mu_t = c(1e-6, 1e6))


    ## find root
    rst        <- uniroot(ft, var_range)$root
    names(rst) <- var_solve

    rst
}

#' Simulate Negative Binomial
#'
#' Simulate negative binomial by poisson-gamma
#'
#' @param n sample size
#' @param mu mean of the negative binomial
#' @param k  dispersion parameter
#' @param t  follow up time. Default value 1.
#'
#' @export
#'
stb_rc_simu_pgamma <- function(n, mu, k, t = 1, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    smp_lambda <- rgamma(n, k, k / mu)
    rst        <- rpois(n, smp_lambda * t)

    if (!is.null(seed))
        set.seed(old_seed)

    data.frame(id     = 1:n,
               lambda = smp_lambda,
               y      = rst)
}


#' Simulate Negative Binomial with Event time
#'
#' Simulate negative binomial by event time
#'
#' @param n sample size
#' @param mu mean of the negative binomial
#' @param k  dispersion parameter
#' @param n_event number of event time
#'
#' @export
#'
stb_rc_simu_events <- function(n, mu, k = 1, n_event = 1, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    f_single <- function(x) {
        lambda  <- rgamma(1, k, k) * mu
        ts      <- rexp(n_event, lambda)
        cum_ts  <- cumsum(ts)
        t_start <- c(0,      cum_ts)
        t_end   <- c(cum_ts, Inf)

        data.frame(id     = x,
                   inx    = seq_len(n_event + 1),
                   lambda = lambda,
                   start  = t_start,
                   end    = t_end,
                   time   = t_end - t_start)
    }

    rst <- lapply(1:n, function(x) f_single(x))

    if (!is.null(seed))
        set.seed(old_seed)

    rbindlist(rst)
}

#' Censor Recurrent Event Data
#'
#'
#' @param dat_rc recurrent events
#' @param t follow up time
#'
#' @export
#'
stb_rc_simu_events_censor <- function(dat_rc, t = 1) {
    dat_rc %>%
        filter(start < t) %>%
        mutate(censor = if_else(end <= t, 0, 1),
               end    = if_else(end <= t, end, t),
               time   = end - start)
}


#' Convert Recurrent Events to NB
#'
#'
#' @param dat_rc recurrent events data frame
#' @param t follow up time
#'
#' @export
#'
stb_rc_simu_events_to_nb <- function(dat_rc, t = 1) {

    dat <- stb_rc_simu_events_censor(dat_rc, t = t)

    ## minus 1 since there is always a censored row
    dat <- dat %>%
        group_by(id) %>%
        summarize(y = n() - 1)

    rst <- dat_rc %>%
        select(id, lambda) %>%
        distinct() %>%
        left_join(dat, by = "id")

    rst
}


#' Calculate covariance based on NB
#'
#' Calculate covariance between recurrent event times
#'
#' @export
#'
stb_rc_nb_cov <- function(mu, k) {
    a       <- k
    b       <- k / mu

    rst_cov <- b^2 / (a - 1)^2 / (a - 2)
    rst_var <- rst_cov + b^2 / (a - 1) / (a - 2)
    rst_cor <- rst_cov / rst_var
    rst_mu  <- b / (a - 1)

    c(cov = rst_cov,
      cor = rst_cor,
      var = rst_var,
      mu  = rst_mu)
}


#' Plot recurrent events
#'
#' Create clinical timeline plot for each patient
#'
#' @param cut_t Cut time point for the analysis
#' @param subset_n Number of patients to be included in the plot
#'
#' @export
#'
stb_rc_plt_fu <- function(dat,
                          t        = 0,
                          subset_n = 0,
                          seed     = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)


    ## censor
    if (t > 0) {
        dat <- stb_rc_simu_events_censor(dat, t = t)
    }

    ## subset
    if (subset_n > 0) {
        id <- dat %>%
            select(id) %>%
            distinct()

        smps <- sample(seq_len(nrow(id)), subset_n)
        id   <- id[smps, ]
        dat  <- id %>%
            left_join(dat, by = "id")
    }


    if (!is.null(seed))
        set.seed(old_seed)

    ## plot
    dat_e1 <- dat %>%
        filter(1 == inx & 0 == censor) %>%
        arrange(time) %>%
        group_by(Arm) %>%
        mutate(x = row_number())

    dat <- dat %>%
        left_join(dat_e1 %>%
                  select(x, Arm, id),
                  by = c("Arm", "id")) %>%
        mutate(is_first = if_else(1 == inx,
                                  "First",
                                  "After First")) %>%
        filter(0 == censor)

    ggplot(data = dat, aes(x = x, y = time)) +
        geom_bar(stat = "identity",
                 aes(fill = is_first)) +
        geom_point(data = dat,
                   aes(x = x, y = end)) +
        coord_flip() +
        theme_bw() +
        theme(legend.position = "none") +
        facet_wrap(~ Arm) +
        labs(y = "Time to Event", x = "Patients")
}
