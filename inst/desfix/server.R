options(shiny.maxRequestSize = 200*1024^2)

library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(shiny)
library(rstan)
library(stringr)
library(parallel)
require(plotly)


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


    output$table_FIX = renderTable({table_FIX()})

    output$density_m = renderPlot({
        density_m() + coord_cartesian(xlim = input$slider_m_zoom)
    })
    output$density_cv = renderPlot({
        density_cv() + coord_cartesian(xlim = input$slider_cv_zoom)
    })
    output$density_y_tilde = renderPlot({
        density_y_tilde() + coord_cartesian(xlim = input$slider_y_tilde_zoom)
    })

    output$pdf_p_prior = renderPlot({pdf_p_prior()})
    output$pdf_p_post = renderPlot({pdf_p_post()})

    output$dose_action = renderTable({dose_action()})

    output$plot_scenario = renderPlot({plot_scenario()})
    output$table_results = renderTable({apply(table_results(),2,as.character)})

    output$plot_scenario_sim = renderPlot({plot_scenario_sim()})
    output$density_y_sim = renderPlot({density_y_sim()})
    output$density_log_y_sim = renderPlot({density_log_y_sim()})
    output$table_results_sim = renderTable({apply(table_results_sim(),2,as.character)})


    ##-----------------slider--------------
    output$slider_m = renderUI({
        return(sliderInput("slider_m",
                           "Threshold of Mean of FIX Functional Activity (%)",
                           min=0,max=200,step=5,value=c(35,70),width="90%"))
    })
    output$slider_cv = renderUI({
        return(sliderInput("slider_cv",
                           "Threshold of Coefficient of Variation of FIX Functional Activity",
                           min=0,max=3,step=0.1,value=c(0.5,1.2),width="90%"))
    })
    output$slider_y_tilde = renderUI({
        return(sliderInput("slider_y_tilde",
                           "Threshold of Predicted FIX Functional Activity (%)",
                           min=0,max=200,step=5,value=c(5,150),width="90%"))
    })

    output$slider_m_zoom = renderUI({
        return(sliderInput("slider_m_zoom", "Zoom in Mean of FIX Functional Activity (%)",
                           min=0,max=round(max(fit_FIX()$post_m)),step=5,
                           value=c(0,round(max(fit_FIX()$post_m))),width="90%"))
    })
    output$slider_cv_zoom = renderUI({
        return(sliderInput("slider_cv_zoom",
                           "Zoom in Coefficient of Variation of FIX Functional Activity",
                           min=0,max=round(max(fit_FIX()$post_cv)),step=0.1,
                           value=c(0,round(max(fit_FIX()$post_cv))),width="90%"))
    })
    output$slider_y_tilde_zoom = renderUI({
        return(sliderInput("slider_y_tilde_zoom",
                           "Zoom in Predicted FIX Functional Activity (%)",
                           min=0,max=round(max(fit_FIX()$post_y_tilde)),step=5,
                           value=c(0,round(max(fit_FIX()$post_y_tilde))),width="90%"))
    })

})
