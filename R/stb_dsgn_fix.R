## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for dose selection fix studies
##
##  DATE:
##      AUGUST, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#' @export
#'
desfix_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Dose selection FIX studies \n\n")
    cat("Design Parameters:\n")
    cat("    ar_dose:    acceptable range of dose\n")
    cat("    TO BE ADDED \n")
}


#' Default design parameter
#'
#'
internal_desfix_dpara <- function() {
    list(n_pt_max             = 6,
         n_pt_min             = 3,
         dlt_thresh_r         = 0.15,
         dlt_thresh_c         = 0.8,
         dlt_prior            = c(0.1, 0.9),
         fix_thresh           = 0.1,
         fix_prior_meanraw    = c(20, 100),
         fix_prior_cv         = c(0.4, 1.2),
         fix_interval_ind     = c(5, 150),
         fix_interval_meanraw = c(35, 70),
         fix_prior_meaninc    = c(0, 50),
         algorithm            = desfix_algorithm_1,
         mean_raw             = c(20, 50, 80, 100),
         cv                   = 0.5,
         dlt_rates            = c(0.01, 0.01, 0.02, 0.02),
         ar_dose              = 2,
         bayes_model          = "independent",
         n_min_borrow         = 3)
}

#' Generate data
#'
#'
#'
desfix_gen_data <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    rst <- desfix_generate_patients_cv(mean_raw   = lst_design$mean_raw,
                                       cv         = lst_design$cv,
                                       dlt_rate   = lst_design$dlt_rates,
                                       n_patients = lst_design$n_pt_max)

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}

#' Simulate patients
#'
#' @param cv Coefficient of variation. If it is scalar, the CV is assumed to be
#'     the same across dose levels. If it is a vector, it has to be the same
#'     length as mean_raw
#'
#' @export
#'
#'
desfix_generate_patients_cv <- function(mean_raw, cv, ...) {

    if (1 == length(cv)) {
        cv <- rep(cv, length(mean_raw))
    } else {
        stopifnot(length(cv) == length(mean_raw))
    }

    mean_log <- log(mean_raw / sqrt(1 + cv^2))
    sd_log   <- sqrt(log(1 + cv^2))

    desfix_generate_patients(mean_log, sd_log, ...)
}

#' Simulate patients
#'
#' the lengths of dose_level, mean_log, sd_log and dlt_rate must be the same
#'
#' right now let n_patients to be the same across dose levels
#'
#' @export
#'
desfix_generate_patients <- function(mean_log, sd_log, dlt_rate, n_patients) {

    stopifnot(length(mean_log) == length(sd_log) &
                  length(mean_log) == length(dlt_rate))

    func_act <- mapply(rlnorm,
                       n       = n_patients,
                       meanlog = mean_log,
                       sdlog   = sd_log)

    dlt      <- mapply(rbinom,
                       n    = n_patients,
                       size = 1,
                       prob = dlt_rate)

    data     <- as.data.frame(cbind(reshape2::melt(func_act)[, 3],
                                    reshape2::melt(dlt)[, 3],
                                    rep(seq_len(length(mean_log)),
                                        each = n_patients)))

    colnames(data) <- c("functional_activity", "dlt", "dose_level")
    data
}


## --------------------------------------------------------------------------
##               INTERNAL FUNCTIONS
## --------------------------------------------------------------------------

#' Get probabilities for intervals
#'
#' @export
#'
get_prob_intervals <- function(post_smps, interval) {
    interval <- c(-Inf, interval, Inf)
    rst      <- sapply(interval,
                       function(x) mean(post_smps <= x))

    sapply(2:length(rst), function(x) rst[x] - rst[x - 1])
}


#' Design point
#'
#' Whether the dose is safe
#'
get_guarded <- function(dta_fix, dta_dlt,
                        fix_interval_ind = c(5, 150),
                        dlt_thresh_r     = 0.15,
                        dlt_thresh_c     = 0.8,
                        dlt_prior        = c(0.1, 0.9)) {


    ## bayesian DLT threshold
    post_a   <- dlt_prior[1] + sum(dta_dlt)
    post_b   <- dlt_prior[2] + sum(1 - dta_dlt)
    p_dlt_gt <- 1 - pbeta(dlt_thresh_r, post_a, post_b)

    ## fix activity and dlt
    rst <- c(all(dta_fix < fix_interval_ind[2]),
             p_dlt_gt < dlt_thresh_c)

    rst
}


#' Dose escalation algorithm
#'
#' 0:stay
#'
#' -12:stop based on predicted and posterior mean FIX given 6 subjects
#' -13:stop based on predicted and posterior mean FIX given 3~5 subjects
#' -14:stop due to both over posterior DLT given >=2 subjects and over observed FIX for any subject
#' -15:stop due to over observed FIX for any subjects
#' -16:stop due to over posterior DLT given >=2 subjects
#' No stop rule for DLT when only 1 subject; Has stop rule for FIX for any subject.
#'
#' 11: escalate due to very low observed FIX of the 1st subjects
#' 12: escalate based on predicted and posterior mean FIX given 6 subjects
#' 13: escalate based on predicted and posterior mean FIX given 3~5 subjects
#'
#' 2:selected based on predicted and posterior mean FIX given 6 subjects
#'
#' @export
#'
desfix_algorithm_1 <- function(cur_data, lst_design, data_hist = NULL, ...) {

    n_pt_max          <- lst_design$n_pt_max
    n_pt_min          <- lst_design$n_pt_min

    dlt_thresh_r      <- lst_design$dlt_thresh_r
    dlt_thresh_c      <- lst_design$dlt_thresh_c
    dlt_prior         <- lst_design$dlt_prior

    fix_thresh        <- lst_design$fix_thresh
    fix_prior_meanraw <- lst_design$fix_prior_meanraw
    fix_prior_cv      <- lst_design$fix_prior_cv
    fix_prior_meaninc <- lst_design$fix_prior_meaninc

    fix_interval_ind       <- lst_design$fix_interval_ind
    fix_interval_meanraw   <- lst_design$fix_interval_meanraw
    bayes_model            <- lst_design$bayes_model

    ## minimum pts to be borrowed from
    n_min_borrow           <- lst_design$n_min_borrow

    dta_fix <- cur_data[, 1]
    dta_dlt <- cur_data[, 2]
    npt      <- nrow(cur_data)
    res      <- rep(NA, npt)

    ## being accelerated?
    if (is.null(data_hist)) {
        is_acc <- 1
    } else {
        tmp_des <- data_hist$decision
        is_acc  <- 0 == length(which(na.omit(tmp_des) != 11))
    }


    if (dta_fix[1] < fix_interval_ind[1] &
        is_acc) {
        res[1] <- 11
        return(res)
    } else if (dta_fix[1] >= fix_interval_ind[2]) { # 1st output of get_guarded()
        res[1] <- -15
        return(res)
    } else {
        res[1] <- 0
    }

    # remove unused simulated data
    if (is.null(data_hist)) {
        data_hist_no_na <- data_hist
    } else {
        data_hist_no_na <- data_hist %>%
            filter(!is.na(decision))
    }

    ## enroll and treat the rest patients
    for (i in 2 : npt) {
        cumu_fix   <- dta_fix[1:i]
        cumu_dlt   <- dta_dlt[1:i]

        ## check dlt and over fix
        is_guarded <- get_guarded(cumu_fix, cumu_dlt,
                                  fix_interval_ind,
                                  dlt_thresh_r,
                                  dlt_thresh_c,
                                  dlt_prior)

        if (all(0 == is_guarded)) {
            res[i] <- -14
            break
        } else if (0 == is_guarded[1]) {
            res[i] <- -15
            break
        } else if (0 == is_guarded[2]) {
            res[i] <- -16
            break
        }

        ## dont check until the first cohort is finished
        if (i < n_pt_min) {
            res[i] <- 0
            next
        }

        ## check posterior fix
        bayes_samples <- desfix_bayes(y            = cumu_fix,
                                      L_m          = fix_prior_meanraw[1],
                                      U_m          = fix_prior_meanraw[2],
                                      L_cv         = fix_prior_cv[1],
                                      U_cv         = fix_prior_cv[2],
                                      L_m_inc      = fix_prior_meaninc[1],
                                      U_m_inc      = fix_prior_meaninc[2],
                                      data_hist    = data_hist_no_na,
                                      n_min_borrow = n_min_borrow,
                                      bayes_model  = bayes_model,
                                      ...)

        post_m      <- bayes_samples$post_m[,ncol(bayes_samples$post_m)]
        pmu         <- get_prob_intervals(post_m, fix_interval_meanraw)
        post_y      <- bayes_samples$post_y_tilde[,ncol(bayes_samples$post_y_tilde)]
        py          <- get_prob_intervals(post_y, fix_interval_ind)

        max_mu      <- which.max(pmu)
        is_over_fix <- py[3] > fix_thresh

        if (n_pt_max == i) {
            if (!is_over_fix & 1 == max_mu) {
                res[i] <- 12
            } else if (!is_over_fix & 2 == max_mu) {
                res[i] <- 2
            } else {
                res[i] <- -12
            }
        } else {
            if (!is_over_fix & 1 == max_mu) {
                res[i] <- 13
                break
            } else if (is_over_fix & 3 == max_mu) {
                res[i] <- -13
                break
            } else {
                res[i] <- 0
            }
        }
    }

    return(res)
}


#' Bayesian log normal models
#'
#' For y, L_m and U_m
#'
#' @export
#'
desfix_bayes = function(y, L_m, U_m, L_cv, U_cv, L_m_inc = NA, U_m_inc = NA,
                        data_hist = NULL,
                        n_min_borrow = 3,
                        bayes_model = c("independent", "same_cv", "monotone"),
                        ...){

    ## check parameters
    if (!is.vector(y))
        stop("y must be a vector.")

    if (length(L_m)  != 1 |
        length(U_m)  != 1 |
        length(L_cv) != 1 |
        length(U_cv) != 1 |
        length(L_m_inc)  != 1 |
        length(U_m_inc)  != 1 )
        stop("L or U must be a single number.")

    bayes_model <- match.arg(bayes_model)

    ## prepare data. Remove data from DL with n < m_min_borrow
    if (is.null(data_hist)) {
        hist_lvs <- data.frame()
    } else {
        hist_lvs <- data_hist %>%
            group_by(dose_level) %>%
            summarize(n = n()) %>%
            filter(n >= n_min_borrow)
    }

    dl <- NULL
    if ("independent" == bayes_model |
        0 == nrow(hist_lvs)) {
        mdl    <- "fix_ind"
        # mdl <- "independent"
        n_dose <- 1
    } else {
        y_temp <- NULL
        for (i in 1:nrow(hist_lvs)) {
            tmp <- data_hist %>%
                filter(dose_level == hist_lvs$dose_level[i]) %>%
                select(functional_activity)

            dl     <- c(dl,     rep(i, nrow(tmp)))
            y_temp <- c(y_temp, tmp$functional_activity)
        }

        n_dose <- nrow(hist_lvs) + 1
        dl     <- c(dl,     rep(n_dose, length(y)))
        y      <- c(y_temp, y)
        mdl    <- switch(bayes_model,
                         "same_cv"  = "fix_samecv",
                         "monotone" = "fix_mono")
        # mdl <- bayes_model
    }

    lst_data = list(N       = length(y),
                    n_dose  = n_dose,
                    DL      = dl,
                    y       = y,
                    L_m     = L_m,
                    U_m     = U_m,
                    L_cv    = L_cv,
                    U_cv    = U_cv,
                    L_m_inc = L_m_inc,
                    U_m_inc = U_m_inc)

    # if (mdl == "independent") {
    #     fit <- stan(file = 'fix_ind.stan', data = lst_data)
    # } else if (mdl == "same_cv"){
    #     fit <- stan(file = 'fix_samecv.stan', data = lst_data)
    # } else if (mdl == "monotone"){
    #     fit <- stan(file = 'fix_mono.stan', data = lst_data)
    # }

    fit <- stb_stan(lst_data, stan_mdl = mdl, ...)

    post_par     = rstan::extract(fit)
    post_m       = data.frame(post_par$m)
    post_cv      = data.frame(post_par$cv)
    post_y_tilde = data.frame(post_par$y_tilde)
    rst <- list(post_m       = post_m,
                post_cv      = post_cv,
                post_y_tilde = post_y_tilde)

    return(rst)
}


#' Single Trial
#'
#'
#' @export
#'
desfix_single_trial <- function(dta_all, lst_design, ...) {

    n_dose  <- max(dta_all$dose_level)
    rst     <- NULL
    for (i in 1:n_dose) {
        cur_data <- dta_all %>%
            dplyr::filter(dose_level == i)

        data_hist <- rst
        cur_rst   <- lst_design$algorithm(cur_data,
                                          lst_design,
                                          data_hist = data_hist,
                                          ...)
        cur_data$decision <- cur_rst
        rst               <- rbind(rst, cur_data)

        ##no escalation
        if (!any(10 < cur_rst, na.rm = TRUE))
            break
    }

    ## remaining patients
    if (i < n_dose) {
        to_add <- dta_all %>%
            dplyr::filter(dose_level > i)
        to_add$decision <- NA

        rst <- rbind(rst, to_add)
    }

    rst
}


#' Summarize results
#'
#' Dose escalation algorithm
#'
#' 0:stay
#'
#' -12:stop based on predicted and posterior mean FIX given 6 subjects
#' -13:stop based on predicted and posterior mean FIX given 3~5 subjects
#' -14:stop due to both over posterior DLT given >=2 subjects and over observed FIX for any subject
#' -15:stop due to over observed FIX for any subjects
#' -16:stop due to over posterior DLT given >=2 subjects
#' No stop rule for DLT when only 1 subject; Has stop rule for FIX for any subject.
#'
#' 11: escalate due to very low observed FIX of the 1st subjects
#' 12: escalate based on predicted and posterior mean FIX given 6 subjects
#' 13: escalate based on predicted and posterior mean FIX given 3~5 subjects
#'
#' 2:selected based on predicted and posterior mean FIX given 6 subjects
#'
#'
#' @export
#'
desfix_summary <- function(results, lst_design, ...) {


    m_true <- lst_design$mean_raw
    n_dose <- length(m_true)

    results$dose_level = factor(results$dose_level,
                                levels = 1 : n_dose)

    n_rep                   = max(results$rep)
    prob_sel                = table(results[results$decision == 2, c("dose_level")]) / n_rep
    prob_stop_after         = table(results[results$decision == -12, c("dose_level")]) / n_rep
    prob_stop_before        = table(results[results$decision == -13, c("dose_level")]) / n_rep
    prob_stop_both_fix_dlt  = table(results[results$decision == -14, c("dose_level")]) / n_rep
    prob_stop_fix           = table(results[results$decision == -15, c("dose_level")]) / n_rep
    prob_stop_dlt           = table(results[results$decision == -16, c("dose_level")]) / n_rep
    prob_es_after_1         = table(results[results$decision == 11, c("dose_level")]) / n_rep
    prob_es_after           = table(results[results$decision == 12, c("dose_level")]) / n_rep
    prob_es_before          = table(results[results$decision == 13, c("dose_level")]) / n_rep

    prob_stop     = sum(results$decision<0)/n_rep
    prob_more_h   = sum(results$decision %in% c(11,12,13) & results$dose_level == n_dose)/n_rep
    ave_n         = apply(table(results[, c("dose_level","rep")]), 1, mean)
    ave_n_lower   = table(results[,"dose_level"],
                          factor(results$functional_activity<lst_design$fix_interval_ind[1],
                                 levels=c(FALSE,TRUE),labels=c(FALSE,TRUE)))[,2]/n_rep
    ave_n_upper   = table(results[,"dose_level"],
                          factor(results$functional_activity>=lst_design$fix_interval_ind[2],
                                 levels=c(FALSE,TRUE),labels=c(FALSE,TRUE)))[,2]/n_rep
    ave_n_total   = mean(apply(table(results[,c("dose_level","rep")]), 2, sum))

    table_results <- data.frame(DLT_rate = lst_design$dlt_rates,
                                m_true   = m_true,
                                prob_sel = as.numeric(round(prob_sel*100,1)),
                                prob_stop_after = as.numeric(round(prob_stop_after*100,1)),
                                prob_stop_before = as.numeric(round(prob_stop_before*100,1)),
                                prob_stop_both_fix_dlt = as.numeric(round(prob_stop_both_fix_dlt*100,1)),
                                prob_stop_fix = as.numeric(round(prob_stop_fix*100,1)),
                                prob_stop_dlt = as.numeric(round(prob_stop_dlt*100,1)),
                                prob_es_after_1 = as.numeric(round(prob_es_after_1*100,1)),
                                prob_es_after = as.numeric(round(prob_es_after*100,1)),
                                prob_es_before = as.numeric(round(prob_es_before*100,1)),
                                ave_n = round(ave_n, 2),
                                ave_n_lower = round(ave_n_lower, 2),
                                ave_n_upper = round(ave_n_upper, 2),
                                ave_n_total=c(rep("",length(m_true)-1),round(ave_n_total,2)),
                                prob_stop = c(rep("",length(m_true)-1),round(prob_stop*100,1)),
                                prob_more_h=c(rep("",length(m_true)-1),round(prob_more_h*100,1)))

    colnames(table_results) = c("DLT Rate",
                                "Mean of FIX (%)",
                                "Prob. Select (%)",
                                paste("Prob. stop based on predicted and posterior mean FIX given",lst_design$n_pt_max,"subjects (%)"),
                                paste("Prob. stop based on predicted and posterior mean FIX given",lst_design$n_pt_min,"~",lst_design$n_pt_max-1,"subjects (%)"),
                                paste("Prob. stop due to both over posterior DLT given >= 2 subjects and over observed FIX for any subject (%)"),
                                paste("Prob. stop due to over observed FIX for any subject (%)"),
                                paste("Prob. stop due to over posterior DLT given >= 2 subjects (%)"),
                                paste("Prob. escalate due to very low observed FIX of the 1st patient (%)"),
                                paste("Prob. escalate based on predicted and posterior mean FIX given",lst_design$n_pt_max,"subjects (%)"),
                                paste("Prob. escalate based on predicted and posterior mean FIX given",lst_design$n_pt_min,"~",lst_design$n_pt_max-1,"subjects (%)"),
                                "Ave. No. of treated subjects",
                                paste("Ave. No. of treated subjects with FIX <",lst_design$fix_interval_ind[1],"%"),
                                paste("Ave. No. of treated subjects with FIX >=",lst_design$fix_interval_ind[2],"%"),
                                "Ave. total No. of treated subjects",
                                "Prob. stop (%)",
                                paste("Prob. select >", n_dose, "(%)"))
    DL_names = rep(NA,n_dose)
    for (i in 1:n_dose){
        DL_names[i] = paste0("DL",i)
    }
    rownames(table_results) = DL_names

    return(t(table_results))
}

## --------------------------------------------------------------------------
##               PRESENTATION FUNCTIONS
## --------------------------------------------------------------------------

#' PLOT DESIGN SIMULATE SCENARIO
#'
#' @export
#'
desfix_plot_scenario <- function(lst_design, ...) {

    dlt_rates            <- lst_design$dlt_rates
    mean_raw             <- lst_design$mean_raw
    fix_interval_meanraw <- lst_design$fix_interval_meanraw
    ar_dose              <- lst_design$ar_dose

    doses     <- seq_len(length(dlt_rates))
    data_plot <- data.frame(Prob = c(dlt_rates*100,
                                     mean_raw,
                                     rep(fix_interval_meanraw[1],
                                         length(dlt_rates)),
                                     rep(fix_interval_meanraw[2],
                                         length(dlt_rates))),
                            DL   = doses,
                            Type = factor(rep(1:4,
                                              each = length(dlt_rates))))

    levels(data_plot$Type) <- c("DLT Rate (%)",
                                "Mean of FIX Activity (%)",
                                "Lower Bound of Target Mean (%)",
                                "Upper Bound of Target Mean (%)")

    plot_scenario <- ggplot(data_plot,
                            aes(x = DL, y = Prob, group = Type)) +
        geom_line(aes(color = Type, linetype = Type), size = 1) +
        geom_point(aes(shape = Type, color = Type), size=4) +
        geom_rect(aes(xmin = min(ar_dose) - 0.2,
                      xmax = max(ar_dose) + 0.2,
                      ymin = 0, ymax = 100,
                      fill = ""),
                  alpha = 0,
                  color = "green3") +
        scale_shape_manual(values = c(19, 19, NA, NA)) +
        scale_color_manual(values = c("red","deepskyblue2","grey","black")) +
        scale_fill_manual(values = "green3",labels="Ideal Dose") +
        ylim(0, 120) +
        ylab("Probability") + xlab("Dose Level") +
        scale_x_continuous(breaks=seq(min(data_plot$DL),max(data_plot$DL),1)) +
        theme_bw() +
        theme(text=element_text(size=20),legend.position="top",
              legend.title=element_blank(),
              legend.key.height= unit(1, 'cm'),
              legend.key.width= unit(2, 'cm')) +
        guides(fill=guide_legend(ncol=1),color=guide_legend(ncol=1),
               shape=guide_legend(ncol=1))

    return(plot_scenario)
}

#' Plot binomial density
#'
#' @export
#'
desfix_beta_binomial_desity = function(x,n,alpha,beta,cut,prior_post){

    if (prior_post=="prior"){
        p_prior_mean = round(alpha/(alpha+beta),3)
        cum_f = round(1-pbeta(cut,alpha,beta),3)
        graph_labels = c(bquote("With Prior Mean of"~p==.(p_prior_mean)),
                         bquote(Pr(p>=.(cut))==.(cum_f)))
        ylab = "Prior Density"
    }

    if (prior_post=="post"){
        p_post_mean = round((alpha+x)/(alpha+beta+n),3)
        cum_f = round(1-pbeta(cut,alpha+x,beta+n-x),3)
        graph_labels = c(bquote("With Posterior Mean of"~p==.(p_post_mean)),
                         bquote(Pr(p>=.(cut)~"| data")==.(cum_f)))
        ylab = "Posterior Density"
    }

    p = seq(0,1,0.001)
    if (prior_post=="prior") f = dbeta(p,alpha,beta)
    if (prior_post=="post") f = dbeta(p,alpha+x,beta+n-x)
    data = data.frame(p=p,f=f,group="1",group2="1")
    data_2 = data.frame(p=p[p>=cut],f=f[p>=cut],group="2",group2="2")
    data = data[data$f!=Inf,]
    data_2 = data_2[data_2$f!=Inf,]

    pdf = ggplot(data,aes(x=p, y=f, fill=group, color=group2)) +
        geom_ribbon(data_2,
                    mapping=aes(ymax=f,fill=group,color=group2),ymin=0,alpha=0.5) +
        geom_line() +
        scale_color_manual(values=c("black","white"),name="",labels=graph_labels) +
        scale_fill_manual(values=c("white","red"),name="",labels=graph_labels) +
        labs(x="DLT Rate (p)",y=ylab) +
        theme_bw() + theme(legend.position="top",text=element_text(size=20))

    return(pdf)
}

#' Plot kernel density
#'
#' @export
#'
desfix_kernel_density = function(x,cut,xlab_name,para_notation,percent,from,to){

    cum_prob_L = round(mean(x<cut[1]),3)
    cum_prob_M = round(mean(x>=cut[1] & x<=cut[2]),3)
    cum_prob_U = round(mean(x>cut[2]),3)

    if (percent=="Y"){
        post_mean = round(mean(x),1)
        graph_labels = c(bquote("With Mean of"~.(para_notation)*"%"==.(post_mean)*"%"),
                         bquote(Pr(.(para_notation)*"%"<.(cut[1])*"%")==.(cum_prob_L)),
                         bquote(Pr(.(cut[1])*"%"<={.(para_notation)*"%"<=.(cut[2])}*"%")==.(cum_prob_M)),
                         bquote(Pr(.(para_notation)*"%">.(cut[2])*"%")==.(cum_prob_U)) )
    } else if (percent=="N"){
        post_mean = round(mean(x),3)
        graph_labels = c(bquote("With Mean of"~.(para_notation)==.(post_mean)),
                         bquote(Pr(.(para_notation)<.(cut[1]))==.(cum_prob_L)),
                         bquote(Pr(.(cut[1])<={.(para_notation)<=.(cut[2])})==.(cum_prob_M)),
                         bquote(Pr(.(para_notation)>.(cut[2]))==.(cum_prob_U)) )
    }

    color = c("black","white","white","white")
    fill = c("white","yellow3","green4","red3")

    if (is.na(from) & is.na(to)) {
        den = density(x)
    } else if (is.na(from) & !is.na(to)){
        den = density(x,to=to)
    } else if (!is.na(from) & is.na(to)){
        den = density(x,from=from)
    } else {
        den = density(x,from=from,to=to)
    }
    data = data.frame(y=den$x,f=den$y,group="1",group2="1")

    data_L = data[data$y<cut[1],]
    data_M = data[data$y>=cut[1] & data$y<=cut[2],]
    data_U = data[data$y>cut[2],]

    ii = 0
    if (dim(data_L)[1]==0) {
        fill=fill[-2-ii];
        color=color[-2-ii];
        graph_labels=graph_labels[-2-ii];
        ii=ii+1
    } else {
        data_L$group="2"; data_L$group2="2"
    }

    if (dim(data_M)[1]==0) {
        fill=fill[-(3-ii)];
        color=color[-3-ii];
        graph_labels=graph_labels[-3-ii];
        ii=ii+1
    } else {
        data_M$group="3"; data_M$group2="3"
    }

    if (dim(data_U)[1]==0) {
        fill=fill[-(4-ii)];
        color=color[-4-ii];
        graph_labels=graph_labels[-4-ii];
        ii=ii+1
    } else {
        data_U$group="4"; data_U$group2="4"
    }

    pdf = ggplot(data,aes(x=y,y=f,fill=group,color=group2)) +
        geom_ribbon(data_L,mapping=aes(ymax=f,fill=group,color=group2),ymin=0) +
        geom_ribbon(data_M,mapping=aes(ymax=f,fill=group,color=group2),ymin=0,
                    show.legend=FALSE) +
        geom_ribbon(data_U,mapping=aes(ymax=f,fill=group,color=group2),ymin=0,
                    show.legend=FALSE) +
        geom_line() +
        scale_color_manual(values=color,name="",labels=graph_labels) +
        scale_fill_manual(values=alpha(fill,0.5),name="",labels=graph_labels) +
        labs(x=xlab_name,y="Density") +
        theme_bw() + theme(legend.position="top",text=element_text(size=20)) +
        guides(fill=guide_legend(nrow=2,byrow=TRUE),color=guide_legend(nrow=2,byrow=TRUE))
    return(pdf)
}



#' Plot multiple kernel density
#'
#' @export
#'
desfix_kernel_density_multiple = function(x,cut,xlab_name,para_notation,logYN,percent,from,to){

    adjust = 4 # make graphs smoother

    color = c("black","white","white","white")
    fill = c("white","yellow3","green4","red3")

    if (logYN=="Y"){
        x$functional_activity = log(x$functional_activity)
    }
    if (logYN!="Y" & logYN!="N"){
        stop("logYN must be Y or N.")
    }

    cum_prob_L = unique(x %>% group_by(dose_level) %>%
                            mutate(p=round(mean(functional_activity<cut[1]),3)) %>%
                            select(dose_level,p))[["p"]]
    cum_prob_M = unique(x %>% group_by(dose_level) %>%
                            mutate(p=round(mean(functional_activity>=cut[1] & functional_activity<=cut[2]),3)) %>%
                            select(dose_level,p))[["p"]]
    cum_prob_U = unique(x %>% group_by(dose_level) %>%
                            mutate(p=round(mean(functional_activity>cut[2]),3)) %>%
                            select(dose_level,p))[["p"]]

    cum_prob_L_label = paste0("[", paste(format(cum_prob_L, nsmall = 3), collapse=", "), "]")
    cum_prob_M_label = paste0("[", paste(format(cum_prob_M, nsmall = 3), collapse=", "), "]")
    cum_prob_U_label = paste0("[", paste(format(cum_prob_U, nsmall = 3), collapse=", "), "]")

    if (percent=="Y"){
        graph_labels = c("Density Curve",
                         bquote(Pr(.(para_notation)*"%"<.(cut[1])*"%")==.(cum_prob_L_label)),
                         bquote(Pr(.(cut[1])*"%"<={.(para_notation)*"%"<=.(cut[2])*"%"})==.(cum_prob_M_label)),
                         bquote(Pr(.(para_notation)*"%">.(cut[2])*"%")==.(cum_prob_U_label)))
    } else if (percent=="N"){
        graph_labels = c("Density Curve",
                         bquote(Pr(.(para_notation)<.(cut[1]))==.(cum_prob_L_label)),
                         bquote(Pr(.(cut[1])<={.(para_notation)<=.(cut[2])})==.(cum_prob_M_label)),
                         bquote(Pr(.(para_notation)>.(cut[2]))==.(cum_prob_U_label)))
    }


    range_DL = range(x$dose_level)
    n_counts = table(x$dose_level)

    den_data = function(x,range_DL,from,to){
        temp_den = sapply(min(range_DL):max(range_DL),function(i){
            if (dim(x[x$dose_level==i,])[1]<2){
                return(data.frame(y=NA,f=NA,dose_level=NA))
            }

            if (is.na(from) & is.na(to)) {
                den_i = density(x[x$dose_level==i,]$functional_activity,adjust=adjust)
            } else if (is.na(from) & !is.na(to)){
                den_i = density(x[x$dose_level==i,]$functional_activity,to=to,adjust=adjust)
            } else if (!is.na(from) & is.na(to)){
                den_i = density(x[x$dose_level==i,]$functional_activity,from=from,adjust=adjust)
            } else {
                den_i = density(x[x$dose_level==i,]$functional_activity,from=from,to=to,adjust=adjust)
            }
            return(data.frame(y=den_i$x,f=den_i$y,dose_level=i))
        })

        den = c()
        for (i in min(range_DL):max(range_DL)){
            den = rbind(den,data.frame(temp_den[,i]))
        }
        return(den)
    }
    data = data.frame(den_data(x,range_DL,from,to),group="1",group2="1")
    data = data[!is.na(data$dose_level),]

    data_L = data[data$y<cut[1],]
    data_M = data[data$y>=cut[1] & data$y<=cut[2],]
    data_U = data[data$y>cut[2],]

    ii = 0
    if (dim(data_L)[1]==0) {
        fill=fill[-2-ii];
        color=color[-2-ii];
        graph_labels=graph_labels[-2-ii];
        ii=ii+1
    } else {
        data_L$group="2"; data_L$group2="2"
    }

    if (dim(data_M)[1]==0) {
        fill=fill[-(3-ii)];
        color=color[-3-ii];
        graph_labels=graph_labels[-3-ii];
        ii=ii+1
    } else {
        data_M$group="3"; data_M$group2="3"
    }

    if (dim(data_U)[1]==0) {
        fill=fill[-(4-ii)];
        color=color[-4-ii];
        graph_labels=graph_labels[-4-ii];
        ii=ii+1
    } else {
        data_U$group="4"; data_U$group2="4"
    }

    labels_DL = c()
    for (i in min(range_DL):max(range_DL)){
        labels_DL = c(labels_DL,paste0("DL ",i," (n = ",n_counts[i],")"))
    }

    data$dose_level = factor(data$dose_level,levels=min(range_DL):max(range_DL),labels=labels_DL)
    data_L$dose_level = factor(data_L$dose_level,levels=min(range_DL):max(range_DL),labels=labels_DL)
    data_M$dose_level = factor(data_M$dose_level,levels=min(range_DL):max(range_DL),labels=labels_DL)
    data_U$dose_level = factor(data_U$dose_level,levels=min(range_DL):max(range_DL),labels=labels_DL)

    pdf = ggplot(data,aes(x=y,y=f,fill=group,color=group2)) +
        facet_wrap(~dose_level, nrow=1) +
        geom_ribbon(data_L,mapping=aes(ymax=f,fill=group,color=group2,group=dose_level),ymin=0) +
        geom_ribbon(data_M,mapping=aes(ymax=f,fill=group,color=group2,group=dose_level),ymin=0,
                    show.legend=FALSE) +
        geom_ribbon(data_U,mapping=aes(ymax=f,fill=group,color=group2,group=dose_level),ymin=0,
                    show.legend=FALSE) +
        geom_line() +
        scale_color_manual(values=color,name="",label=graph_labels) +
        scale_fill_manual(values=alpha(fill,0.5),name="",label=graph_labels) +
        labs(x=xlab_name,y="Density") +
        theme_bw() + theme(legend.position="top",text=element_text(size=20)) +
        guides(fill=guide_legend(nrow=4,byrow=TRUE),color=guide_legend(nrow=4,byrow=TRUE))
    return(pdf)
}


#' Plot multiple density for FIX
#'
#' @export
#'
lognormal_density_multiple = function(cut,true_m,true_cv,xlab_name,para_notation,logYN,percent){

    color = c("black","white","white","white")
    fill = c("white","yellow3","green4","red3")

    mu = log(true_m / sqrt(1 + true_cv^2))
    if (length(true_cv)==1) {
        sigma = rep(sqrt(log(1 + true_cv^2)),length(true_m))
    } else {
        sigma = sqrt(log(1 + true_cv^2))
    }

    if (logYN=="Y"){
        cum_prob_L = round(pnorm(cut[1],mu,sigma),3)
        cum_prob_M = round(pnorm(cut[2],mu,sigma)-pnorm(cut[1],mu,sigma),3)
        cum_prob_U = round(1-pnorm(cut[2],mu,sigma),3)
        y = seq(min(mu-3*sigma),max(mu+3*sigma),0.05)
    } else if (logYN=="N"){
        cum_prob_L = round(plnorm(cut[1],mu,sigma),3)
        cum_prob_M = round(plnorm(cut[2],mu,sigma)-plnorm(cut[1],mu,sigma),3)
        cum_prob_U = round(1-plnorm(cut[2],mu,sigma),3)
        y = seq(0,exp(max(mu+3*sigma)),1)
    }

    cum_prob_L_label = paste0("[", paste(format(cum_prob_L, nsmall = 3), collapse=", "), "]")
    cum_prob_M_label = paste0("[", paste(format(cum_prob_M, nsmall = 3), collapse=", "), "]")
    cum_prob_U_label = paste0("[", paste(format(cum_prob_U, nsmall = 3), collapse=", "), "]")

    if (percent=="Y"){
        graph_labels = c("Density Curve",
                         bquote(Pr(.(para_notation)*"%"<.(cut[1])*"%")==.(cum_prob_L_label)),
                         bquote(Pr(.(cut[1])*"%"<={.(para_notation)*"%"<=.(cut[2])*"%"})==.(cum_prob_M_label)),
                         bquote(Pr(.(para_notation)*"%">.(cut[2])*"%")==.(cum_prob_U_label)))
    } else if (percent=="N"){
        graph_labels = c("Density Curve",
                         bquote(Pr(.(para_notation)<.(cut[1]))==.(cum_prob_L_label)),
                         bquote(Pr(.(cut[1])<={.(para_notation)<=.(cut[2])})==.(cum_prob_M_label)),
                         bquote(Pr(.(para_notation)>.(cut[2]))==.(cum_prob_U_label)))
    }

    data = c()
    labels_DL = c()
    for (i in 1:length(true_m)){
        if (logYN=="Y") {
            f_y = dnorm(y,mu[i],sigma[i])
        } else if (logYN=="N") {
            f_y = dlnorm(y,mu[i],sigma[i])
        }
        data = rbind(data,data.frame(y=y,f_y=f_y,dose_level=i,group="1",group2="1"))
        labels_DL = c(labels_DL,paste0("DL ",i))
    }
    data$dose_level = factor(data$dose_level,levels=1:length(true_m),labels=labels_DL)

    data_L = data[data$y<cut[1],]
    data_M = data[data$y>=cut[1] & data$y<=cut[2],]
    data_U = data[data$y>cut[2],]

    ii = 0
    if (dim(data_L)[1]==0) {
        fill=fill[-2-ii];
        color=color[-2-ii];
        graph_labels=graph_labels[-2-ii];
        ii=ii+1
    } else {
        data_L$group="2"; data_L$group2="2"
    }

    if (dim(data_M)[1]==0) {
        fill=fill[-(3-ii)];
        color=color[-3-ii];
        graph_labels=graph_labels[-3-ii];
        ii=ii+1
    } else {
        data_M$group="3"; data_M$group2="3"
    }

    if (dim(data_U)[1]==0) {
        fill=fill[-(4-ii)];
        color=color[-4-ii];
        graph_labels=graph_labels[-4-ii];
        ii=ii+1
    } else {
        data_U$group="4"; data_U$group2="4"
    }

    pdf = ggplot(data,aes(x=y,y=f_y,fill=group,color=group2)) +
        facet_wrap(~dose_level, nrow=1) +
        geom_ribbon(data_L,mapping=aes(ymax=f_y,fill=group,color=group2,group=dose_level),ymin=0) +
        geom_ribbon(data_M,mapping=aes(ymax=f_y,fill=group,color=group2,group=dose_level),ymin=0,
                    show.legend=FALSE) +
        geom_ribbon(data_U,mapping=aes(ymax=f_y,fill=group,color=group2,group=dose_level),ymin=0,
                    show.legend=FALSE) +
        geom_line() +
        scale_color_manual(values=color,name="",label=graph_labels) +
        scale_fill_manual(values=alpha(fill,0.5),name="",label=graph_labels) +
        labs(x=xlab_name,y="Density") +
        theme_bw() + theme(legend.position="top",text=element_text(size=20)) +
        guides(fill=guide_legend(nrow=4,byrow=TRUE),color=guide_legend(nrow=4,byrow=TRUE))
    return(pdf)
}
