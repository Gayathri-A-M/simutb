## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX GENERIC FUNCTIONS
##
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##                  class stb_design and class stb_dist
## -----------------------------------------------------------------------------

setGeneric("stb_describe",
           function(x, ...) standardGeneric("stb_describe"))

setGeneric("stb_get_para",
           function(x) standardGeneric("stb_get_para"))

setGeneric("stb_para<-",
           function(x, value) standardGeneric("stb_para<-"))

setGeneric("stb_set_default_para",
           function(x, ...) standardGeneric("stb_set_default_para"))

setGeneric("stb_generate_data",
           function(x, ...) standardGeneric("stb_generate_data"))

setGeneric("stb_analyze_data",
           function(x, data, ...) standardGeneric("stb_analyze_data"))

setGeneric("stb_plot_data",
           function(x, data, ...) standardGeneric("stb_plot_data"))

setGeneric("stb_create_trial",
           function(x, ...) standardGeneric("stb_create_trial"))

setGeneric("stb_create_simustudy",
           function(x,
                    n_rep  = 5,
                    n_core = 5,
                    seed   = NULL, ...) standardGeneric("stb_create_simustudy"))

setGeneric("stb_simu_gen_raw",
           function(x, lst, ...) standardGeneric("stb_simu_gen_raw"))

setGeneric("stb_simu_gen_summary",
           function(x, lst, ...) standardGeneric("stb_simu_gen_summary"))

setGeneric("stb_simu_gen_key",
           function(x, lst, ...) standardGeneric("stb_simu_gen_key"))


## -----------------------------------------------------------------------------
##                        class stb_trial
## -----------------------------------------------------------------------------

setGeneric("stb_get_trial_data",
           function(x) standardGeneric("stb_get_trial_data"))

setGeneric("stb_get_trial_result",
           function(x) standardGeneric("stb_get_trial_result"))

setGeneric("stb_get_trial_seed",
           function(x) standardGeneric("stb_get_trial_seed"))

setGeneric("stb_get_trial_design",
           function(x) standardGeneric("stb_get_trial_design"))

setGeneric("stb_trial_plot",
           function(x, ...) standardGeneric("stb_trial_plot"))


## -----------------------------------------------------------------------------
##                        class stb_simustudy
## -----------------------------------------------------------------------------

setGeneric("stb_get_simu_design",
           function(x) standardGeneric("stb_get_simu_design"))

setGeneric("stb_get_simu_raw",
           function(x) standardGeneric("stb_get_simu_raw"))

setGeneric("stb_get_simu_summary",
           function(x) standardGeneric("stb_get_simu_summary"))

setGeneric("stb_get_simu_key",
           function(x) standardGeneric("stb_get_simu_key"))

setGeneric("stb_get_simu_nrep",
           function(x) standardGeneric("stb_get_simu_nrep"))

setGeneric("stb_get_simu_seed",
           function(x) standardGeneric("stb_get_simu_seed"))


## -----------------------------------------------------------------------------
##                        class stb_dist
## -----------------------------------------------------------------------------

setGeneric("stb_set_para_by_qs",
           function(x, target_qs, ...) standardGeneric("stb_set_para_by_qs"))

setGeneric("stb_set_para_by_qs_ind",
           function(x, target_qs, ...)
               standardGeneric("stb_set_para_by_qs_ind"))

setGeneric("stb_dist_get_qs",
           function(x, qs) standardGeneric("stb_dist_get_qs"))


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       HELPER FUNCTIONS
##
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

#' Create a distribution
#'
#'
#'
#' @export
#'
stb_dist_by_qs <- function(type = c("normal",
                                    "lognormal",
                                    "beta"),
                           target_qs = c("0.5" = 0.7,
                                         "0.8" = 0.9)) {

    type <- match.arg(type)
    rst  <- switch(type,
                   normal    = new("STB_DIST_NORMAL"),
                   lognormal = new("STB_DIST_LOGNORMAL"),
                   beta      = new("STB_DIST_BETA"),
                   new("STB_DIST"))

    stb_set_para_by_qs(rst, target_qs = target_qs)
}

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
                                       "bayes_2arm")) {

    type <- match.arg(type)
    rst  <- switch(type,
                   surv_strat = new("STB_DESIGN_STRAT_SURV"),
                   surv_or    = new("STB_DESIGN_SURV_OR"),
                   surv_join  = new("STB_DESIGN_SURV_JOIN"),
                   surv_biom  = new("STB_DESIGN_SURV_BIOM"),
                   bayes_1arm = new("STB_DESIGN_BAYES_1ARM"),
                   bayes_2arm = new("STB_DESIGN_BAYES_2ARM"),
                   new("STB_DESIGN"))

    rst
}

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
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##                        overall class stb_design
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN",
         slots     = list(design_para  = "list",
                          design_valid = "numeric"),
         prototype = prototype(design_para  = list(),
                               design_valid = 1))

setMethod("initialize",
          "STB_DESIGN",
          function(.Object, ...) {
              .Object@design_para <- stb_set_default_para(.Object)
              .Object
          })

setMethod("stb_set_default_para",
          "STB_DESIGN",
          function(x) list())

setMethod("stb_get_para",
          "STB_DESIGN",
          function(x) x@design_para)

setMethod("stb_para<-",
          "STB_DESIGN",
          function(x, value) {
              x@design_para <- tl_merge_lists(value,
                                              x@design_para)
              x
          })

setMethod("stb_describe",
          "STB_DESIGN",
          function(x, ...) NULL)

setMethod("stb_generate_data",
          "STB_DESIGN",
          function(x, ...) NULL)

setMethod("stb_analyze_data",
          "STB_DESIGN",
          function(x, data, ...) list())

setMethod("stb_create_trial",
          "STB_DESIGN",
          function(x, seed = NULL, ...) {

              if (!is.null(seed))
                  old_seed <- set.seed(seed)

              data   <- stb_generate_data(x, ...)
              result <- stb_analyze_data(x, data)

              if (!is.null(seed))
                  set.seed(old_seed)

              new("STB_TRIAL",
                  design = x,
                  data   = data,
                  result = result,
                  seed   = seed)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN",
          function(x, lst, ...) list())

setMethod("stb_simu_gen_summary",
          "STB_DESIGN",
          function(x, lst, ...) list())

setMethod("stb_simu_gen_key",
          "STB_DESIGN",
          function(x, lst, ...) list())


setMethod("stb_create_simustudy",
          "STB_DESIGN",
          function(x,
                   n_rep  = 5,
                   n_core = 5,
                   seed   = NULL, ...) {

              if (0 == x@design_valid) {
                  cat("Invalid design. \n")
                  return(NULL)
              }

              ## seed
              if (!is.null(seed))
                  old_seed <- set.seed(seed)

              ## all random seeds
              all_seeds <- ceiling(abs(rnorm(n_rep) * 100000))

              ## replications
              rst <-
                  parallel::mclapply(
                                seq_len(n_rep),
                                function(k) {
                                    if (0 == k %% 5)
                                        print(k)

                                    cur_trial <-
                                        stb_create_trial(
                                            x,
                                            seed = all_seeds[k],
                                            ...)

                                    cur_trial@result
                                },
                                mc.cores = n_core)

              ## summarize
              rst_raw     <- stb_simu_gen_raw(x, rst)
              rst_summary <- stb_simu_gen_summary(x, rst_raw)
              rst_key     <- stb_simu_gen_key(x, rst_summary)

              ## seed
              if (!is.null(seed))
                  set.seed(old_seed)

              new("STB_SIMU_STUDY",
                  design      = x,
                  n_rep       = n_rep,
                  n_core      = n_core,
                  rst_raw     = rst_raw,
                  rst_summary = rst_summary,
                  rst_key     = rst_key,
                  seed        = seed)
          })


## -----------------------------------------------------------------------------
##                        overall class stb_trial
## -----------------------------------------------------------------------------
#'
#' @export
#'
setClass("STB_TRIAL",
         slots = list(design = "STB_DESIGN",
                      data   = "data.frame",
                      result = "list",
                      seed   = "numeric"),
         prototype = prototype(seed   = NULL,
                               result = NULL))

setMethod("stb_get_trial_data",   "STB_TRIAL", function(x) x@data)
setMethod("stb_get_trial_result", "STB_TRIAL", function(x) x@result)
setMethod("stb_get_trial_design", "STB_TRIAL", function(x) x@design)
setMethod("stb_get_trial_seed",   "STB_TRIAL", function(x) x@seed)

setMethod("stb_trial_plot",
          "STB_TRIAL",
          function(x, ...) {
              stb_plot_data(x@design, x@data, ...)
          })

## -----------------------------------------------------------------------------
##                        overall class stb_simustudy
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_SIMU_STUDY",
         slots = list(design      = "STB_DESIGN",
                      n_rep       = "numeric",
                      n_core      = "numeric",
                      rst_raw     = "list",
                      rst_summary = "list",
                      rst_key     = "list",
                      seed        = "numeric"))


setMethod("stb_get_simu_design",   "STB_SIMU_STUDY", function(x) x@design)
setMethod("stb_get_simu_raw",      "STB_SIMU_STUDY", function(x) x@rst_raw)
setMethod("stb_get_simu_summary",  "STB_SIMU_STUDY", function(x) x@rst_summary)
setMethod("stb_get_simu_key",      "STB_SIMU_STUDY", function(x) x@rst_key)
setMethod("stb_get_simu_nrep",     "STB_SIMU_STUDY", function(x) x@n_rep)
setMethod("stb_get_simu_seed",     "STB_SIMU_STUDY", function(x) x@seed)


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
          function(x, data, ...) {
              callNextMethod()
              rst_interim <- survjoin_ana_data(data,
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
          function(x, data, ...) {
              callNextMethod()
              survbiom_ana_data(data, x@design_para)
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
          function(x, data, ...) {
              callNextMethod()
              rst_post <- bayes_ana_data(
                  data,
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

setClass("STB_DESIGN_BAYES_1ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_1ARM",
          function(x) {
              internal_bayes1arm_dpara()
          })

setClass("STB_DESIGN_BAYES_2ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_2ARM",
          function(x) {
              internal_bayes2arm_dpara()
          })


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX CLASSES OF DISTRIBUTIONS
##
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##                        overall class stb_dist
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DIST",
         slots     = list(dist_para  = "list"),
         prototype = prototype(dist_para  = list()))

setMethod("stb_get_para",
          "STB_DIST",
          function(x) x@dist_para)

setMethod("stb_para<-",
          "STB_DIST",
          function(x, value) {
              x@dist_para <- tl_merge_lists(value,
                                            x@dist_para)
              x
          })

setMethod("stb_set_para_by_qs",
          "STB_DIST",
          function(x,
                   target_qs = c("0.1" = 0.5,
                                 "0.7" = 0.9),
                   ...) {

              para  <- stb_set_para_by_qs_ind(x,
                                              target_qs = target_qs,
                                              ...)
              para$target_qs <- target_qs
              x@dist_para    <- para
              x
          })

setClass("STB_DIST_NORMAL",
         contains  = "STB_DIST",
         slots     = list(dist_para  = "list"),
         prototype = prototype(dist_para  = list(mean  = 0,
                                                 sd    = 1)))
setMethod("stb_set_para_by_qs_ind",
          "STB_DIST_NORMAL",
          function(x, target_qs, ...)
              stb_dist_normal_para_by_qs(target_qs, ...))

setMethod("stb_dist_get_qs",
          "STB_DIST_NORMAL",
          function(x, qs) {
              rst <- qnorm(qs,
                           mean = x@dist_para$mean,
                           sd   = x@dist_para$sd)
              names(rst) <- qs
              rst
          })


setClass("STB_DIST_BETA",
         contains  = "STB_DIST",
         slots     = list(dist_para  = "list"),
         prototype = prototype(dist_para  = list(shape1 = 1,
                                                 shape2 = 1)))
setMethod("stb_set_para_by_qs_ind",
          "STB_DIST_BETA",
          function(x, target_qs, ...)
              stb_dist_beta_para_by_qs(target_qs, ...))

setMethod("stb_dist_get_qs",
          "STB_DIST_BETA",
          function(x, qs) {
              rst <- qbeta(qs,
                           shape1 = x@dist_para$shape1,
                           shape2 = x@dist_para$shape2)
              names(rst) <- qs
              rst})

setClass("STB_DIST_LOGNORMAL",
         contains  = "STB_DIST",
         slots     = list(dist_para  = "list"),
         prototype = prototype(dist_para  = list(mean = 1,
                                                 sd   = 1)))

setMethod("stb_set_para_by_qs_ind",
          "STB_DIST_LOGNORMAL",
          function(x, target_qs, ...)
              stb_dist_lognormal_para_by_qs(target_qs, ...))

setMethod("stb_dist_get_qs",
          "STB_DIST_LOGNORMAL",
          function(x, qs) {
              rst <- qnorm(qs,
                           mean = x@dist_para$norm_mean,
                           sd   = x@dist_para$norm_sd)

              rst        <- exp(rst)
              names(rst) <- qs
              rst})
