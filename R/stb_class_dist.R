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
##       HELPER FUNCTIONS
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

#'
#' @export
#'
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


#'
#' @export
#'
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

#'
#' @export
#'
setClass("STB_DIST_LOGNORMAL",
         contains  = "STB_DIST",
         slots     = list(dist_para  = "list"),
         prototype = prototype(dist_para  = list(mean = 1,
                                                 sd   = 1)))

#'
#' @export
#'
setMethod("stb_set_para_by_qs_ind",
          "STB_DIST_LOGNORMAL",
          function(x, target_qs, ...)
              stb_dist_lognormal_para_by_qs(target_qs, ...))

#'
#' @export
#'
setMethod("stb_dist_get_qs",
          "STB_DIST_LOGNORMAL",
          function(x, qs) {
              rst <- qnorm(qs,
                           mean = x@dist_para$norm_mean,
                           sd   = x@dist_para$norm_sd)

              rst        <- exp(rst)
              names(rst) <- qs
              rst})
