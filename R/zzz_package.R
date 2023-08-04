#' The 'simutb' package.
#'
#' @docType package
#' @name    simutb-package
#' @aliases simutb
#' @useDynLib simutb, .registration = TRUE
#'
#' @import methods
#' @import stats
#' @import ggplot2
#'
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines par plot points text arrows grid
#'     rect
#' @importFrom parallel detectCores
#' @importFrom utils as.roman
#' @importFrom dplyr %>% group_by_ group_by summarize mutate count mutate_if
#'     rename filter select arrange ungroup n distinct left_join if_else rowwise
#' @importFrom tidyr gather
#' @importFrom data.table rbindlist
#' @importFrom survival Surv survfit coxph survdiff
#' @importFrom rpact getDesignGroupSequential
#' @importFrom mvtnorm rmvnorm pmvnorm dmvnorm
#' @importFrom survminer ggsurvplot ggsurvplot_facet
#' @importFrom mmrm mmrm df_1d
#' @importFrom brms brm as_draws_matrix
#' @importFrom scales percent
#'
NULL
