## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX CLASSES OF DESIGNS
##
##
##   DESIGNS:
##      1. STB_DESIGN_STRAT_SURV
##      2. STB_DESIGN_SURV_OR
##      3. STB_DESIGN_SURV_JOIN
##      4. STB_DESIGN_SURV_BIOM
##      5. STB_DESIGN_BAYES_1ARM
##      6. STB_DESIGN_BAYES_2ARM
##      7. STB_DESIGN_RCURRENT
##      8. STB_DESIGN_RMEASURE
##      9. STB_DESIGN_MSMA_SURV
##     10: STB_DESIGN_DOSE_FIX
##     11: STB_DESIGN_COVID
##     12: STB_DESIGN_RCURRENT_ADAPT
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##       HELPER FUNCTIONS
## -----------------------------------------------------------------------------

#' Create a Study Design
#'
#'
#'
#' @export
#'
stb_create_design <- function(type = c("surv_strat",
                                       "surv_or",
                                       "surv_join",
                                       "surv_biom",
                                       "bayes_1arm",
                                       "bayes_2arm",
                                       "rcurrent",
                                       "rmeasure",
                                       "msma_surv",
                                       "dose_fix",
                                       "covid",
                                       "rcurrent_adapt")) {

    type <- match.arg(type)
    rst  <- switch(type,
                   surv_strat     = new("STB_DESIGN_STRAT_SURV"),
                   surv_or        = new("STB_DESIGN_SURV_OR"),
                   surv_join      = new("STB_DESIGN_SURV_JOIN"),
                   surv_biom      = new("STB_DESIGN_SURV_BIOM"),
                   bayes_1arm     = new("STB_DESIGN_BAYES_1ARM"),
                   bayes_2arm     = new("STB_DESIGN_BAYES_2ARM"),
                   rcurrent       = new("STB_DESIGN_RCURRENT"),
                   rmeasure       = new("STB_DESIGN_RMEASURE"),
                   msma_surv      = new("STB_DESIGN_MSMA_SURV"),
                   dose_fix       = new("STB_DESIGN_DOSE_FIX"),
                   covid          = new("STB_DESIGN_COVID"),
                   rcurrent_adapt = new("STB_DESIGN_RCURRENT_ADAPT"),
                   new("STB_DESIGN"))

    rst
}


## -----------------------------------------------------------------------------
##                        stratified survival
## -----------------------------------------------------------------------------
setClass("STB_DESIGN_STRAT_SURV",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_STRAT_SURV",
          function(x, ...) {
              callNextMethod()
              stratsurv_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_STRAT_SURV",
          function(x) {
              stratsurv_default_para()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_STRAT_SURV",
          function(x, ...) {
              if (0 == x@design_valid) {
                  cat("Invalid design \n")
                  return(NULL)
              }

              strasurv_gen_data(x@design_para, ...)
          })


## -----------------------------------------------------------------------------
##                        oncology trial with OR
## -----------------------------------------------------------------------------
setClass("STB_DESIGN_SURV_OR",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_SURV_OR",
          function(x, ...) {
              callNextMethod()
              survor_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_SURV_OR",
          function(x) {
              survor_default_para()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_SURV_OR",
          function(x, ...) {
              if (0 == x@design_valid) {
                  cat("Invalid design \n")
                  return(NULL)
              }

              do.call(survor_gen_data,
                      c(x@design_para, ...))
          })


## -----------------------------------------------------------------------------
##                        surival with joint PFS and OS
## -----------------------------------------------------------------------------
setClass("STB_DESIGN_SURV_JOIN",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_SURV_JOIN",
          function(x, ...) {
              callNextMethod()
              survjoin_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_SURV_JOIN",
          function(x) {
              survjoin_default_para()
          })

setMethod("stb_para<-",
          "STB_DESIGN_SURV_JOIN",
          function(x, value) {
              x             <- callNextMethod(x, value)
              lst_para      <- do.call(internal_survjoin_dpara,
                                       x@design_para)
              x@design_para <- lst_para
              x
          })

setMethod("stb_generate_data",
          "STB_DESIGN_SURV_JOIN",
          function(x, ...) {
              if (0 == x@design_valid) {
                  cat("Invalid design \n")
                  return(NULL)
              }

              do.call(survjoin_gen_data, c(x@design_para, ...))
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_SURV_JOIN",
          function(x, data_ana, ...) {
              callNextMethod()
              rst_interim <- survjoin_ana_data(data_ana[[1]],
                                               x@design_para$target_primary,
                                               x@design_para$info_frac,
                                               x@design_para$primary)

              rst <- survjoin_hiertest(rst_interim,
                                       x@design_para$pval_bound)

              rst
          })

setMethod("stb_plot_data",
          "STB_DESIGN_SURV_JOIN",
          function(x, data, ...) {
              survjoin_plot_data(data, ...)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_SURV_JOIN",
          function(x, lst, ...) {
              rst <- list()
              for (i in seq_len(length(lst))) {
                  cur_lst     <- lst[[i]]
                  cur_lst$rep <- i
                  rst[[i]]    <- cur_lst
              }

              list(rbindlist(rst))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_SURV_JOIN",
          function(x, lst, ...) {
              rst <- survjoin_simu_summary(
                  lst[[1]],
                  colnames(x@design_para$pval_bound))
              rst
          })

setMethod("stb_simu_gen_key",
          "STB_DESIGN_SURV_JOIN",
          function(x, lst, ...) {
              rst <- survjoin_simu_key(lst,
                                       x@design_para)
              list(rst_key = rst)
          })


## -----------------------------------------------------------------------------
##                        surival with arm dropping
## -----------------------------------------------------------------------------
setClass("STB_DESIGN_SURV_BIOM",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_SURV_BIOM",
          function(x, ...) {
              callNextMethod()
              survbiom_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_SURV_BIOM",
          function(x) {
              survbiom_default_para()
          })

setMethod("stb_para<-",
          "STB_DESIGN_SURV_BIOM",
          function(x, value) {
              x             <- callNextMethod(x, value)
              lst_para      <- do.call(internal_survbiom_dpara,
                                       x@design_para)
              x@design_para <- lst_para
              x
          })

setMethod("stb_generate_data",
          "STB_DESIGN_SURV_BIOM",
          function(x, ...) {
              if (0 == x@design_valid) {
                  cat("Invalid design \n")
                  return(NULL)
              }

              survbiom_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_SURV_BIOM",
          function(x, data_ana, ...) {
              callNextMethod()
              survbiom_ana_data(data_ana[[1]], x@design_para)
          })

setMethod("stb_plot_data",
          "STB_DESIGN_SURV_BIOM",
          function(x, data, ...) {
              survbiom_plot_data(data, ...)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_SURV_BIOM",
          function(x, lst, ...) {

              n_reps      <- length(lst)
              rst_biom    <- list()
              rst_interim <- list()
              for (i in seq_len(length(lst))) {
                  rst_biom[[i]]    <- lst[[i]]$biom %>% mutate(rep = i)
                  rst_interim[[i]] <- lst[[i]]$rst  %>% mutate(rep = i)
              }

              list(n_reps      = n_reps,
                   rst_biom    = rbindlist(rst_biom),
                   rst_interim = rbindlist(rst_interim))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_SURV_BIOM",
          function(x, lst, ...) {
              rst <- survbiom_simu_summary(lst)
              rst
          })

setMethod("stb_simu_gen_key",
          "STB_DESIGN_SURV_BIOM",
          function(x, lst, ...) {
              rst <- survbiom_simu_key(lst,
                                       x@design_para)
              list(rst_key = rst)
          })


## -----------------------------------------------------------------------------
##                        bayesian 1- and 2-arm design
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_BAYES",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_BAYES",
          function(x, ...) {
              callNextMethod()
              bayes_describe(x, ...)
          })

setMethod("stb_generate_data",
          "STB_DESIGN_BAYES",
          function(x, ...) {
              bayes_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_BAYES",
          function(x, data_ana, ...) {
              callNextMethod()
              rst_post <- bayes_ana_data(
                  data_ana[[1]],
                  x@design_para$prior_by_arm,
                  n_post = x@design_para$n_post,
                  x      = x@design_para$x_post,
                  ...)

              rst_decision <- bayes_ana_decision(
                  rst_post,
                  decision_ref    = x@design_para$decision_ref,
                  decision_h0     = x@design_para$decision_h0,
                  decision_gl     = x@design_para$decision_gl,
                  decision_thresh = x@design_para$decision_thresh)

              list(rst_post    = rst_post,
                   rst_diff    = rst_decision$rst_diff,
                   rst_success = rst_decision$rst_success)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_BAYES",
          function(x, lst, ...) {
              rst <- list()
              for (i in seq_len(length(lst))) {
                  cur_lst     <- data.frame(lst[[i]]$rst_success)
                  cur_lst$rep <- i
                  rst[[i]]    <- cur_lst
              }

              list(rbindlist(rst))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_BAYES",
          function(x, lst, ...) {
              rst <- lst[[1]] %>%
                  group_by(arm) %>%
                  summarize(success = mean(success))
              list(rst)
          })

#'
#' @export
#'
setClass("STB_DESIGN_BAYES_1ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_1ARM",
          function(x) {
              internal_bayes1arm_dpara()
          })

#'
#' @export
#'
setClass("STB_DESIGN_BAYES_2ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_2ARM",
          function(x) {
              internal_bayes2arm_dpara()
          })


## -----------------------------------------------------------------------------
##                        recurrent event
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_RCURRENT",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_RCURRENT",
          function(x, ...) {
              callNextMethod()
              rcurrent_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_RCURRENT",
          function(x) {
              internal_rcurrent_dpara()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_RCURRENT",
          function(x, ...) {
              rcurrent_gen_data(x@design_para, ...)
          })

#'
#' @export
#'
setMethod("stb_create_analysis_set",
          "STB_DESIGN_RCURRENT",
          function(x, data,
                   type          = c("min_fu", "fix_fu"),
                   fu_days       = 12 * 7,
                   min_fu_days   = 12 * 7,
                   pt_proportion = 1,
                   ...) {


              if (is.null(data))
                  return(NULL)

              type     <- match.arg(type)
              dat_full <- switch(
                  type,
                  min_fu = rcurrent_day_eos_1(data,
                                              min_fu_days   = min_fu_days,
                                              pt_proportion = pt_proportion),

                  fix_fu = rcurrent_day_eos_2(data,
                                              fu_days       = fu_days,
                                              pt_proportion = pt_proportion)
              )

              dat_full <- rcurrent_censor(dat_full)
              dat_nb   <- rcurrent_get_nb(dat_full)

              list(data    = dat_full,
                   data_nb = dat_nb)
          })


## -----------------------------------------------------------------------------
##                        repeated measure
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_RMEASURE",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_RMEASURE",
          function(x, ...) {
              callNextMethod()
              rmeasure_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_RMEASURE",
          function(x) {
              internal_rmeasure_dpara()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_RMEASURE",
          function(x, ...) {
              rmeasure_gen_data(x@design_para, ...)
          })

setMethod("stb_create_analysis_set",
          "STB_DESIGN_RMEASURE",
          function(x, data, par_interim = NULL, ...) {

              if (is.null(data))
                  return(NULL)

              if (is.null(par_interim))
                  par_interim <- x@design_para$par_interim

              rst <- list()
              for (i in seq_len(nrow(par_interim))) {
                  rst[[i]] <- rmeasure_create_ana_data(data,
                                                       par_interim[i, 1],
                                                       par_interim[i, 2])
              }

              rst
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_RMEASURE",
          function(x, data_ana, par_analysis = NULL, ...) {

              if (is.null(data_ana))
                  return(NULL)

              if (is.null(par_analysis))
                  par_analysis <- x@design_para$par_analysis

              rst <- list()
              for (i in seq_len(length(data_ana))) {
                  cur_d    <- data_ana[[i]]
                  cur_mmrm <- rmeasure_ana_mmrm(
                      cur_d,
                      endpoint_visit = par_analysis$endpoint_visit)
                  cur_ttest <- rmeasure_ana_ttest(
                      cur_d,
                      endpoint_visit = par_analysis$endpoint_visit)
                  cur_bayes <- rmeasure_ana_bayes(
                      cur_d,
                      n_visit = ncol(x@design_para$mu_by_arm),
                      ...)

                  rst[[i]] <- list(mmrm  = cur_mmrm,
                                   ttest = cur_ttest,
                                   bayes = cur_bayes)
              }

              rst
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_RMEASURE",
          function(x, lst, ...) {
              rst <- rmeasure_simu_summary(lst, ...)
              rst
          })

## -----------------------------------------------------------------------------
##                        msma with survival outcome
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_MSMA_SURV",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_MSMA_SURV",
          function(x, ...) {
              callNextMethod()
              msma_surv_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_MSMA_SURV",
          function(x) {
              internal_msma_surv_dpara()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_MSMA_SURV",
          function(x, ...) {
              msma_surv_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_MSMA_SURV",
          function(x, data_ana, par_analysis = NULL, ...) {

              if (is.null(par_analysis))
                  par_analysis <- x@design_para$par_analysis

              list(msma_surv_ana_logrank(data_ana[[1]]))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_MSMA_SURV",
          function(x, lst, ...) {
              rst <- msma_surv_simu_summary(lst, ...)
              rst
          })


## -----------------------------------------------------------------------------
##                        Dose escaltion for FIX study
## -----------------------------------------------------------------------------
#'
#' @export
#'
setClass("STB_DESIGN_DOSE_FIX",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              callNextMethod()
              desfix_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_DOSE_FIX",
          function(x) {
              internal_desfix_dpara()
          })

setMethod("stb_plot_design",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              desfix_plot_scenario(x@design_para, ...)
          })

setMethod("stb_generate_data",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              desfix_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_DOSE_FIX",
          function(x, data_ana, ...) {
              rst <- desfix_single_trial(data_ana[[1]],
                                         lst_design = x@design_para,
                                         ...)

              list(rst)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_DOSE_FIX",
          function(x, lst, ...) {
              n_reps  <- length(lst)
              rst     <- list()
              for (i in seq_len(length(lst))) {
                  rst[[i]]  <- lst[[i]][[1]] %>% mutate(rep = i)
              }

              list(n_reps = n_reps,
                   rst    = rbindlist(rst))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_DOSE_FIX",
          function(x, lst, ...) {
              rst <- desfix_summary(
                  lst$rst, x@design_para, ...)
              list(rst)
          })


## -----------------------------------------------------------------------------
##                        COVID STUDIES
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_COVID",
         contains = "STB_DESIGN_MSMA_SURV")

setMethod("stb_describe",
          "STB_DESIGN_COVID",
          function(x, ...) {
              covid_describe()
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_COVID",
          function(x) {
              lst <- callNextMethod()
              lst <- internal_covid_dpara(lst)
              lst
          })

setMethod("stb_create_analysis_set",
          "STB_DESIGN_COVID",
          function(x, data, par_interim = NULL, ...) {

              if (is.null(data))
                  return(NULL)

              if (is.null(par_interim))
                  par_interim <- x@design_para$par_interim

              rst <- list()
              for (i in seq_len(length(par_interim$ana_fraction))) {
                  rst[[i]] <- stb_tl_interim_data(
                      data,
                      total     = par_interim$target_primary,
                      info_frac = par_interim$ana_fraction[i],
                      event     = "obs",
                      arms      = par_interim$target_arms)
              }

              names(rst) <- par_interim$ana_fraction

              rst
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_COVID",
          function(x, data_ana, par_analysis = NULL, ...) {

              if (is.null(par_analysis))
                  par_analysis <- x@design_para$par_analysis

              rst <- NULL
              for (i in seq_len(length(data_ana))) {

                  if (!is.data.frame(data_ana[[i]]))
                      next

                  info_frac <- names(data_ana)[i]
                  cur_rst   <- covid_analysis(data_ana[[i]],
                                              info_frac,
                                              x@design_para)
                  rst <- rbind(rst, cur_rst)
              }

              list(rst)
          })


setMethod("stb_simu_gen_summary",
          "STB_DESIGN_COVID",
          function(x, lst, ...) {
              rst <- covid_simu_summary(lst, x@design_para$par_analysis, ...)
              rst
          })


## -----------------------------------------------------------------------------
##                   recurrent event with adaptive sample size
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_RCURRENT_ADAPT",
         contains = "STB_DESIGN_RCURRENT")

setMethod("stb_generate_data",
          "STB_DESIGN_RCURRENT_ADAPT",
          function(x, ...) {
              rcurrent_gen_data(x@design_para, ...)
          })

#'
#' @export
#'
setMethod("stb_create_analysis_set",
          "STB_DESIGN_RCURRENT_ADAPT",
          function(x, data, ...) {

              if (is.null(data))
                  return(NULL)

              data_interim <- rcurrent_day_eos_adapt_1(data,
                                                       x@design_para$n_stage1,
                                                       x@design_para$fix_fu)

              data_interim_nb <- rcurrent_get_nb(data_interim)

              ## TO BE ADDED
              ## INPUT: data_interim_nb, alpha, and power
              ## OUTPUT: n_stage2, target_event

              data_final <- rcurrent_day_eos_adapt_2(
                  data,
                  n_stage1 = x@design_para$n_stage1,
                  n_stage2,
                  target_event,
                  rcur_info = x@design_para$rcur_info,
                  fix_fu    = x@design_para$fix_fu)

              data_final_nb <- rcurrent_get_nb(data_final)

              list(data            = data,
                   data_interim    = data_interim,
                   data_interim_nb = data_interim_nb,
                   data_final      = data_final,
                   data_final_nb   = data_final)
          })


setMethod("stb_analyze_data",
          "STB_DESIGN_RCURRENT_ADAPT",
          function(x, data_ana) {

              dat_final    <- data_ana$data_final
              dat_final_nb <- data_ana$data_final_nb

              rst <- stb_tl_rc_reg(dat_final_nb)
              ## sample size and duration
              rst$study_n   <- length(unique(dat_final$sid))
              rst$study_dur <- max(dat_final$date_eos) - max(dat_final$bos)

              list(rst)
          })
