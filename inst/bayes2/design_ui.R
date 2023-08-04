
##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------

tab_main <- function() {

    div(style = "padding-top: 30px;",

        ## Sidebar with a slider input for number of bins
        sidebarLayout(
            sidebarPanel(

                h4("Number of Patients in Treatment Group:"),
                numericInput('trt_n_min', 'Minimum Number of Patients', 40),
                numericInput('trt_n_max', 'Maximum Number of Patients', 120),
                numericInput('trt_n_step', 'Step of Number of Patients', 20),

                h4("Randomization Ratio between Treatment and Control:"),
                numericInput('trt_r','Treatment Group:', 1),
                numericInput('cont_r', 'Control Group:', 1),

                h4("Assumptions:"),
                numericInput('delta', 'Assumed Proportion Difference', 0.15),
                numericInput('trt_cap', 'Treatment Group Rate Cap', 0.95),
                numericInput('simN', 'Simulation Number', 2000),
                selectInput('prior_dist',
                            'Prior Distribution of Placebo Response Rate:',
                            c('Beta Distribution', 'Uniform Distribution')),
                conditionalPanel('input.prior_dist=="Beta Distribution"',
                                 numericInput('prior.alpha',
                                              'Alpha of Prior Distribution',
                                              34.7),
                                 numericInput('prior.beta',
                                              'Beta of Prior Distribution',
                                              13.8)),

                h4("Criteria:"),
                numericInput('alpha', 'Confidence Level alpha:', 0.95),

                ## succ_cri_word <- reactive(paste0('Success Criteria (lower
                ## limit of the one-sided ', input.alpha*100,'% credible
                ## interval for the group difference (treatment - conntrol) is
                ## greater than)')),

                ## numericInput('success_c', succ_cri_word, -0.05),
                numericInput('success_c',
                             'Success Criteria (lower limit of the one-sided 95% credible interval for the group difference (treatment - conntrol) is greater than )',
                             -0.05),
                numericInput('super_c',
                             'Superiority Criteria (lower limit of the one-sided 95% credible interval for the group difference (treatment - conntrol) is greater than )',
                             0)
            ),

            ## Show a plot of the generated distribution
            mainPanel(
                tableOutput("simu_res"),
                br(),
                plotOutput("plot_success"),
                br(),
                plotOutput("plot_superiority")
            )
        )
        )
}


##-------------------------------------------------------------
##           DATA FUNCTIONS
##-------------------------------------------------------------

sim_dd <- reactive({
    trt_n  <- seq(input$trt_n_min,
                  input$trt_n_max,
                  input$trt_n_step)

    cont_n <- round(trt_n * input$cont_r / input$trt_r)

    dt <- lapply(1 : length(trt_n),
                 function(yy) {
        simdata <- lapply(1:input$simN, function(xx) {

            ## Set up per protocol assumption
            if (input$prior_dist == "Beta Distribution") {
                ## Control group response rate: random number from historical
                ## data of Beta distribution
                cont_p <- rbeta(1,
                                input$prior.alpha,
                                input$prior.beta)
            } else if (input$prior_dist == "Uniform Distribution") {
                cont_p <- runif(1, 0, 1)
            }

            ## Aflibercept response rate
            trt_p <- cont_p + input$delta
            ## trt_p > 0.95   # check if treatment group rate > 0.95
            trt_p <- ifelse(trt_p > input$trt_cap,
                            input$trt_cap,
                            trt_p)   ## Give a cap, e.g. 0.95

            ## Create Treatment and Control sample data: 0 - non-responder, 1 -
            ## responder
            cont_dt <- rbinom(cont_n[yy], 1, cont_p)
            trt_dt  <- rbinom(trt_n[yy],  1, trt_p)

            ## Treatment and Control samples' response rates
            trt_phat  <- sum(trt_dt)  / trt_n[yy]
            cont_phat <- sum(cont_dt) / cont_n[yy]

            ## ----------------------------------------------------------------
            ## Sample data Comparison

            ## It should be more straightforward to get samples of treatment and
            ## control response rates from their posterior beta distributions.
            ## Then you can get LCL numerically

            ## Why 5? It should not be fixed in the code and it does not have to
            ## be same as the data generation prior.
            diff     <- trt_phat - cont_phat     # Sample treatment difference
            ## Sample pooled response rate
            pooled_p <- (sum(trt_dt) + sum(cont_dt)) / (trt_n[yy] + cont_n[yy])
            se_CrI   <- sqrt(pooled_p * (1 - pooled_p) *
                             ((1 / (trt_n[yy] + 5)) + (1 / (cont_n[yy] + 5))))
            ## ----------------------------------------------------------------

            ## 95% one-sided credible lower limit
            LCL <- diff - qnorm(input$alpha) * se_CrI

            ## Success criterion: lower limit of the one-sided 95% credible
            ## interval for the treatment difference (treatment – control) is
            ## greater than -5%

            ## 1 - success, 0 - fail
            success <- ifelse(LCL > input$success_c, 1, 0)

            ## superiority criterion: the lower limit of the 95% one-sided
            ## credible interval (treatment – control) > 0

            ## 1 - superiority, 0 - non-superiority
            superiority <- ifelse(LCL > input$super_c, 1, 0)

            ## Summary
            data.frame(diff_phat   = diff,
                       LCL         = LCL,
                       success     = success,
                       superiority = superiority)

        })

        dd          <- data.frame(do.call(rbind, simdata))
        success     <- scales::percent(sum(dd$success)     / input$simN, 0.1)
        superiority <- scales::percent(sum(dd$superiority) / input$simN, 0.1)

        data.frame(trt_n       = trt_n[yy],
                   cont_n      = cont_n[yy],
                   success     = success,
                   superiority = superiority)

    })

    data.frame(do.call(rbind, dt))
})
