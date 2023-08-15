options(shiny.maxRequestSize = 200*1024^2)
require(plotly)
require(simutb)

shinyServer(function(input, output, session) {

    source("design_ui.R", local = TRUE)
    userLog <- reactiveValues()

    ##--------------------------------------
    ##---------main-------------------------
    ##--------------------------------------
    output$mainpage <- renderUI({
        tab_main()
    })

    ##--------------------------------------
    ##---------exit-------------------------
    ##--------------------------------------
    observeEvent(input$close, {
        stopApp()})

    ##--------------------------------------
    ##---------plot-------------------------
    ##--------------------------------------

    output$simu_res <- renderTable({
        sim_dd()
    })

    output$plot_success <- renderPlot({
        ggplot(sim_dd(),
               aes(trt_n, success, group = 1)) +
            geom_line(color = "green") +
            geom_point()
    })

    output$plot_superiority <- renderPlot({
        ggplot(sim_dd(),
               aes(trt_n, superiority, group = 1)) +
            geom_line(color = "blue") +
            geom_point()
    })

})
