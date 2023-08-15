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
         fix_prior_meanraw    = c(0, 300),
         fix_prior_cv         = c(0, 2),
         fix_interval_ind     = c(5, 150),
         fix_interval_meanraw = c(35, 70),
         algorithm            = desfix_algorithm_1,
         mean_raw             = c(20, 50, 80, 100),
         cv                   = 0.5,
         dlt_rates            = c(0.01, 0.01, 0.02, 0.02),
         ar_dose              = 2)
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
desfix_generate_patients <- function(mean_log   = c(2.88, 3.8, 4.27, 4.49),
                                     sd_log     = rep(0.47, 4),
                                     dlt_rate   = c(0.01, 0.01, 0.02, 0.02),
                                     n_patients = 6) {

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

    colnames(data) <- c("functional_activity", "DLT", "dose_level")
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
    post_a   <- dlt_prior[1] + sum(dta_dlt)     + 1
    post_b   <- dlt_prior[1] + sum(1 - dta_dlt) + 1
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
#' -12:stop after 6 patients
#' -13:stop before 6 patients
#' -14:stop for both DLT and Over FIX
#' -15:stop for Over FIX
#' -16:stop for DLT
#'
#' 11: escalate after the 1st patient
#' 12: escalate after 6 patients
#' 13: escalate before 6 patients
#'
#' 2:selected
#'
#' @export
#'
desfix_algorithm_1 <- function(cur_data, lst_design, ...) {

    n_pt_max          <- lst_design$n_pt_max
    n_pt_min          <- lst_design$n_pt_min

    dlt_thresh_r      <- lst_design$dlt_thresh_r
    dlt_thresh_c      <- lst_design$dlt_thresh_c
    dlt_prior         <- lst_design$dlt_prior

    fix_thresh        <- lst_design$fix_thresh
    fix_prior_meanraw <- lst_design$fix_prior_meanraw
    fix_prior_cv      <- lst_design$fix_prior_cv

    fix_interval_ind       <- lst_design$fix_interval_ind
    fix_interval_meanraw   <- lst_design$fix_interval_meanraw


    dta_fix <- cur_data[, 1]
    dta_dlt <- cur_data[, 2]
    npt      <- nrow(cur_data)
    res      <- rep(NA, npt)

    ## check the first patient
    if (dta_fix[1] < fix_interval_ind[1]) {
        res[1] <- 11
        return(res)
    } else {
        res[1] <- 0
    }

    ## enroll and treat the rest patients
    for (i in 2 : npt) {

        ## dont check until the first cohort is finished
        if (i < n_pt_min) {
            res[i] <- 0
            next
        }

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

        ## check posterior fix
        bayes_samples <- desfix_bayes(y    = cumu_fix,
                                      L_m  = fix_prior_meanraw[1],
                                      U_m  = fix_prior_meanraw[2],
                                      L_cv = fix_prior_cv[1],
                                      U_cv = fix_prior_cv[2],
                                      ...)

        pmu  <- get_prob_intervals(bayes_samples$post_m,
                                   fix_interval_meanraw)
        py   <- get_prob_intervals(bayes_samples$post_y_tilde,
                                   fix_interval_ind)

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

        cur_rst  <- lst_design$algorithm(cur_data, lst_design, ...)
        rst      <- c(rst, cur_rst)

        if (!any(10 < cur_rst, na.rm = TRUE))
            break
    }

    to_add <- nrow(dta_all) - length(rst)
    if (to_add > 0)
        rst <- c(rst, rep(NA, to_add))

    dta_all$decision <- rst
    dta_all
}

#' Bayesian log normal models
#'
#'
#' @export
#'
desfix_bayes <- function(y, L_m, U_m, L_cv, U_cv, ...){

    if (!is.vector(y))
        stop("y must be a vector.")

    if (length(L_m)  != 1 |
        length(U_m)  != 1 |
        length(L_cv) != 1 |
        length(U_cv) != 1)
        stop("L or U must be a single number.")

    lst_data <- list(N    = length(y),
                     y    = y   / 100,
                     L_m  = L_m / 100,
                     U_m  = U_m / 100,
                     L_cv = L_cv,
                     U_cv = U_cv)

    fit <- stb_stan(lst_data,
                    stan_mdl = "logn",
                    ...)

    post_par     = rstan::extract(fit)
    post_m       = post_par$m * 100
    post_cv      = post_par$cv
    post_y_tilde = post_par$y_tilde * 100

    rst <- list(post_m       = post_m,
                post_cv      = post_cv,
                post_y_tilde = post_y_tilde)

    return(rst)
}

#' Summarize results
#'
#' @export
#'
desfix_summary <- function(results, lst_design, ...) {

    m_true <- lst_design$mean_raw
    n_dose <- length(m_true)

    results$dose_level = factor(results$dose_level,
                                levels = 1 : n_dose)

    n_rep         = max(results$rep)
    prob_sel      = table(results[results$decision == 2, c("dose_level")]) / n_rep
    prob_stop     = sum(results$decision < 0)/n_rep
    prob_more_h   = sum(results$decision == 1 & results$dose_level == n_dose)/n_rep
    ave_n         = apply(table(results[, c("dose_level","rep")]), 1, mean)
    ave_n_total   = mean(apply(table(results[,c("dose_level","rep")]), 2, sum))

    table_results <- data.frame(Dose     = n_dose,
                                DLT_rate = lst_design$dlt_rates,
                                m_true   = m_true,
                                prob_sel = as.numeric(round(prob_sel*100,1)),
                                ave_n    = round(ave_n, 2),
                                ave_n_total=c(rep("",length(m_true)-1),round(ave_n_total,2)),
                                prob_stop = c(rep("",length(m_true)-1),round(prob_stop*100,1)),
                                prob_more_h=c(rep("",length(m_true)-1),round(prob_more_h*100,1)))

    colnames(table_results) = c("Dose","DLT Rate (%)","Mean of FIX F. A. (%)",
                                "Prob. Sel. (%)","Ave. No. of Subj. per Cohort",
                                "Ave. Total No. of Subj.","Prob. Stop (%)","Prob. Sel.>4 (%)")

    table_results
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
    data_plot <- data.frame(Prob = c(dlt_rates,
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
        ylim(0, 100) +
        ylab("Probability") + xlab("Dose Level") +
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

#' Plot binomial density
#'
#' @export
#'
desfix_kernel_density = function(x,cut,xlab_name,para_notation,percent,from,to){

    cum_prob_L = round(mean(x<cut[1]),3)
    cum_prob_M = round(mean(x>=cut[1] & x<=cut[2]),3)
    cum_prob_U = round(mean(x>cut[2]),3)

    if (percent=="N"){
        post_mean = round(mean(x),3)
        graph_labels = c(bquote("With Mean of"~.(para_notation)==.(post_mean)),
                         bquote(P(.(para_notation)<.(cut[1]))==.(cum_prob_L)),
                         bquote(P(.(cut[1])<={.(para_notation)<=.(cut[2])})==.(cum_prob_M)),
                         bquote(P(.(para_notation)>.(cut[2]))==.(cum_prob_U)) )
    }
    if (percent=="Y"){
        post_mean = round(mean(x),1)
        graph_labels = c(bquote("With Mean of"~.(para_notation)*"%"==.(post_mean)*"%"),
                         bquote(P(.(para_notation)*"%"<.(cut[1])*"%")==.(cum_prob_L)),
                         bquote(P(.(cut[1])*"%"<={.(para_notation)*"%"<=.(cut[2])}*"%")==.(cum_prob_M)),
                         bquote(P(.(para_notation)*"%">.(cut[2])*"%")==.(cum_prob_U)) )
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
        ii=ii+1}
    else {data_L$group="2"; data_L$group2="2"}

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
        data_U$group="4";
        data_U$group2="4"
    }

    data_combine = rbind(data,data_L,data_M,data_U)
    data_combine$group = factor(data_combine$group)
    data_combine$group2 = factor(data_combine$group2)

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
