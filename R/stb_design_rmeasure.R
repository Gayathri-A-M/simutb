## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for repeated measurement studies
##
##  DATE:
##      APRIL, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#' @export
#'
rmeasure_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Two-arm study with repeated measurements \n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:    total sample size (default 100)\n")
    cat("    ratio_by_arm:   randomization proportion for each arm \n")
    cat("                    (default c(0.5, 0.5))\n")
    cat("    mu_by_arm:      mean vector for each arm\n")
    cat("    cov_by_arm:     sigma matrix for each arm\n")
    cat("    visit_day:      visit time in days\n")
    cat("    alpha:          alpha level (default 0.025)\n")
    cat("    par_interim:    list of parameters for creating interim \n")
    cat("                    and final analysis dataset\n")
    cat("    par_analysis:   list of analysis parameters\n")
    cat("    par_enroll:     list of enrollment parameters \n")
}


#' Default design parameter for single arm
#'
#'
internal_rmeasure_dpara <- function() {
    list(sample_size       = 100,
         ratio_by_arm      = c(0.5, 0.5),
         mu_by_arm         = rbind(c(1, 1, 1),
                                   c(1, 1, 1)),
         cov_by_arm        = rbind(c(1,   0.1, 0.1),
                                   c(0.1,   1, 0.1),
                                   c(0.1, 0.1,   1),
                                   c(1,   0.1, 0.1),
                                   c(0.1,   1, 0.1),
                                   c(0.1, 0.1,   1)),
         visit_day         = ceiling(c(6, 9, 12) * 30.25),
         par_enroll        = list(type       = "by_duration",
                                  pt_dur_mth = 16),
         par_interim       = cbind(ana_visit    = c(2 , 3),
                                   ana_fraction = c(0.5, 1)),
         par_analysis      = list(ana_formua = as.formula(
                                      y ~ arm * visit_no +
                                          us(visit_no | sid)),
                                  and_contrast   = NULL,
                                  endpoint_visit = 3),
         alpha             = 0.025
         )
}

#' Generate data
#'
#'
#'
rmeasure_gen_data <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## n by arm
    n_by_arm <- tl_draw_arm_size(sample_size  = lst_design$sample_size,
                                 ratio_by_arm = lst_design$ratio_by_arm)

    ## enrollment
    dat_enroll <- do.call(stb_tl_simu_enroll_arms,
                          list(n_by_arm   = c(n_by_arm),
                               par_enroll = lst_design$par_enroll))

    ## outcome
    ## annual to day hr
    mu_by_arm   <- lst_design$mu_by_arm
    cov_by_arm  <- lst_design$cov_by_arm
    visit_day   <- lst_design$visit_day

    n_visit     <- ncol(mu_by_arm)
    rst         <- NULL
    for (i in seq_len(length(n_by_arm))) {
        cur_y <- stb_tl_simu_rm(n_by_arm[i],
                                mu_by_arm[i, ],
                                cov_by_arm[(i - 1) * n_visit + 1:n_visit, ])


        dim(cur_y) <- c(length(cur_y), 1)
        cur_rst    <- data.frame(arm       = i - 1,
                                 sid       = rep(1:n_by_arm[i], n_visit),
                                 visit_no  = rep(1:n_visit,
                                                 each = n_by_arm[i]),
                                 visit_day = rep(visit_day,
                                                 each = n_by_arm[i]),
                                 y         = cur_y)

        rst <- rbind(rst, cur_rst)
    }

    ## merget data
    rst <- dat_enroll %>%
        left_join(rst, by = c("arm", "sid")) %>%
        mutate(day_visit = day_enroll + visit_day,
               sid       = sid + arm * lst_design$sample_size)


    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}


#' Generate analysis data set
#'
#' Generate analysis dataset for minimum fu days
#'
rmeasure_create_ana_data <- function(data_full, visit_no, fraction) {

    vno      <- visit_no
    dat      <- data_full %>%
        filter(visit_no == vno) %>%
        select(arm, sid, day_visit) %>%
        arrange(day_visit)

    day_eos <- dat[ceiling(nrow(dat) * fraction), "day_visit"]

    data_full %>%
        mutate(day_eos = day_eos) %>%
        filter(day_visit < day_eos)
}


#' Generate analysis data set
#'
#' Generate analysis dataset for minimum fu days
#'
rmeasure_ana_mmrm <- function(data_ana,
                              endpoint_visit = 1,
                              formula        = NULL,
                              contrast       = NULL) {

    dat <- data_ana %>%
        mutate(sid      = factor(sid),
               visit_no = factor(visit_no))

    if (is.null(formula))
        formula <- as.formula(y ~ arm * visit_no + us(visit_no | sid))

    fit <- mmrm::mmrm(formula = formula, data = dat)

    if (is.null(contrast)) {
        contrast <- numeric(length(fit$beta_est))
        names(contrast) <- names(fit$beta_est)
        contrast["arm"] <- 1

        if (endpoint_visit > 1)
            contrast[paste("arm:visit_no", endpoint_visit, sep = "")] <- 1
    }

    ## return
    list(mmrm_fit  = fit,
         mmrm_test = mmrm::df_1d(fit, contrast))
}

rmeasure_ana_ttest <- function(data_ana, endpoint_visit = 1) {

    dat <- data_ana %>%
        filter(visit_no == endpoint_visit)

    t.test(y ~ arm, data = dat)
}

#' Get posterior samples of treatment effect at each visit
#'
#'
rmeasure_ana_bayes <- function(data_ana,
                               n_visit = 3,
                               formula = NULL,
                               ...) {

    dat <- data_ana %>%
        mutate(sid      = factor(sid),
               visit_no = factor(visit_no))

    if (is.null(formula))
        formula <- as.formula(y ~ arm * visit_no + unstr(visit_no, sid))

    fit       <- brms::brm(formula = formula, data = dat, ...)
    post_smps <- brms::as_draws_matrix(fit)

    rst <- post_smps[, "b_arm"]
    for (i in 2:n_visit) {
        cur_rst <- post_smps[, "b_arm"] +
            post_smps[, paste("b_arm:visit_no", i, sep = "")]
        rst     <- cbind(rst, cur_rst)
    }

    rst
}


#' Summarize individual simulation study
rmeasure_simu_ind <- function(lst) {
    rst <- NULL
    for (i in 1:length(lst)) {

        cur_rst <- NULL
        ## MMRM
        cur_mmrm <- lst[[i]]$mmrm
        cur_rst  <- c(i,
                      cur_mmrm$mmrm_test$est,
                      cur_mmrm$mmrm_test$se,
                      cur_mmrm$mmrm_test$p_val)

        ## T Test
        cur_ttest <- lst[[i]]$ttest
        cur_rst   <- c(cur_rst, cur_ttest$p.value)

        ## Bayes
        cur_bayes  <- lst[[i]]$bayes
        n_visit    <- ncol(cur_bayes)
        cur_rst    <- c(cur_rst,
                        apply(cur_bayes, 2, function(x) mean(x > 0)))

        ## bind
        rst <- rbind(rst, cur_rst)
    }

    colnames(rst) <- c("ana_no",
                       "mmrm_est", "mmrm_se", "mmrm_pval",
                       "ttest_pval",
                       paste("p0_visit", 1:n_visit, sep = ""))

    rst
}

#' Summarize simulation study results
#'
#'
rmeasure_simu_summary <- function(lst, ...) {
    rst <- NULL
    for (i in 1:length(lst)) {
        cur_rst     <- rmeasure_simu_ind(lst[[i]], ...)
        cur_rst     <- data.frame(cur_rst)
        cur_rst$rep <- i
        rst         <- rbind(rst, cur_rst)
    }

    list(rst)
}
