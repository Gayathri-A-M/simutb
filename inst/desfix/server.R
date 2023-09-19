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
rstan_options(auto_write = TRUE) 


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

    output$plot_scenario_SR = renderPlot({plot_scenario_SR()})
    output$density_y_true_SR = renderPlot({
      density_y_true_SR() + coord_cartesian(xlim = input$slider_y_true_SR_zoom, 
                                            ylim = c(0,input$slider_f_true_SR_zoom))
    })
    output$density_log_y_true_SR = renderPlot({
      density_log_y_true_SR()} + coord_cartesian(xlim = input$slider_log_y_true_SR_zoom, 
                                                 ylim = c(0,input$slider_f_log_true_SR_zoom))
    )
    output$density_y_sim_SR = renderPlot({
      density_y_sim_SR() + coord_cartesian(xlim = input$slider_y_sim_SR_zoom, 
                                           ylim = c(0,input$slider_f_sim_SR_zoom))
    })
    output$density_log_y_sim_SR = renderPlot({
      density_log_y_sim_SR() + coord_cartesian(xlim = input$slider_log_y_sim_SR_zoom, 
                                               ylim = c(0,input$slider_f_log_sim_SR_zoom))
    })
    output$table_results_SR = renderTable({table_results_SR()},
                                          rownames = TRUE)
    
    output$plot_scenario_CS = renderPlot({plot_scenario_CS()})
    output$density_y_true_CS = renderPlot({
      density_y_true_CS() + coord_cartesian(xlim = input$slider_y_true_CS_zoom,
                                            ylim = c(0,input$slider_f_true_CS_zoom))
    })
    output$density_log_y_true_CS = renderPlot({
      density_log_y_true_CS() + coord_cartesian(xlim = input$slider_log_y_true_CS_zoom,
                                                ylim = c(0,input$slider_f_log_true_CS_zoom))
    })
    output$density_y_sim_CS = renderPlot({
      density_y_sim_CS() + coord_cartesian(xlim = input$slider_y_sim_CS_zoom,
                                           ylim = c(0,input$slider_f_sim_CS_zoom))
    })
    output$density_log_y_sim_CS = renderPlot({
      density_log_y_sim_CS() + coord_cartesian(xlim = input$slider_log_y_sim_CS_zoom,
                                               ylim = c(0,input$slider_f_log_sim_CS_zoom))
    })
    output$table_results_CS = renderTable({table_results_CS()},
                                           rownames = TRUE)
    

})
