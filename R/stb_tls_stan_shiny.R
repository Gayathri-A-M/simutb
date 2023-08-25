## -----------------------------------------------------------------------------
##
##                 STAN
##
## -----------------------------------------------------------------------------

#' Call STAN models
#'
#' Call STAN models. Called by \code{psrwe_powerp}.
#'
#' @param lst_data List of study data to be passed to STAN
#' @param stan_mdl STAN model name
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#' @return Result from STAN sampling
#'
#' @export
#'
stb_stan <- function(lst_data,
                     stan_mdl = c("fix_ind", "fix_samecv", "fix_mono"),
                     chains = 4, iter = 2000, warmup = 1000, cores = 4,
                     control = list(adapt_delta = 0.95), ...) {

    stan_mdl <- match.arg(stan_mdl)
    stan_rst <- rstan::sampling(stanmodels[[stan_mdl]],
                                data    = lst_data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                cores   = cores,
                                control = control,
                                ...)

    stan_rst
}


## -----------------------------------------------------------------------------
##
##                 SHINY
##
## -----------------------------------------------------------------------------

#' Run Web-Based application
#'
#' Call Shiny to run \code{statidea} as a web-based application.
#'
#' @details
#'
#' A web browser will be brought up for users to access the GUI
#'
#'
#' @export
#'
stb_shiny <- function(appname = "bayes2", pkgname = "simutb") {

    req_pkgs        <- c("shiny", "shinythemes", "DT",
                         "knitr", "rmarkdown", "pander")

    chk_uninstalled <- sapply(req_pkgs,
                              function(x) {
                                  !requireNamespace(x,
                                                    quietly = TRUE)
                              })

    chk_inx         <- which(chk_uninstalled)

    if (0 < length(chk_inx)) {

        msg <- paste("For the Shiny app to work, please install ",
                     ifelse(1 < length(chk_inx), "packages ", "package "),
                     paste(req_pkgs[chk_inx], collapse = ", "),
                     " by \n install.packages(",
                     paste(paste("'",
                                 req_pkgs[chk_inx],
                                 "'",
                                 sep = ""), collapse = ", "),
                     ") \n  ",
                     sep = "")

        stop(msg, call. = FALSE)
    }

    app_dir <- system.file(appname, package = "simutb")
    if (app_dir == "") {
        stop("Could not find Shiny directory. Try re-installing `simutb`.",
             call. = FALSE)
    }

    shiny::runApp(app_dir, display.mode = "normal")
}
