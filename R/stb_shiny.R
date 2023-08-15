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
