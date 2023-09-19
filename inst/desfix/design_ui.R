
##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------

##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_design_parameters(),
                tab_conduct_study(),
                tab_simu(),
                tab_simulation_study(),
                tab_other_settings()
                )
}

tab_design_parameters <- function(){
  tabPanel("Design Parameters",
           fluidRow(
             column(6,radioButtons("bayes_model", "Bayesian Model",
                                   c("Independent DLs"="independent","Same CV of FIX"="same_cv", 
                                     "Same CV and Monotone Increasing FIX"="monotone"),inline=TRUE)),
             conditionalPanel(
               condition="input.bayes_model=='same_cv' | input.bayes_model=='monotone'",
               column(3,numericInput("n_min_borrow","Minimum No. of subjects per DL to be borrowed from",3)),
             ),
           ),
           
           fluidRow(
             column(3,numericInput("alpha","Parameter alpha of Beta Prior for DLT Rate",0.1)),
             column(3,numericInput("beta","Parameter beta of Beta Prior for DLT Rate",0.9)),
           ),
           fluidRow(
             column(3,sliderInput("LU_m", "Bounds of Uniform Prior for Mean (m%)",
                                  min=0,max=300,step=5,value=c(20,100))),
             column(3,sliderInput("LU_cv", "Bounds of Uniform Prior for Coefficient of Variation (CV)",
                                  min=0,max=3,step=0.1,value=c(0.4,1.2))),
             conditionalPanel(
               condition="input.bayes_model=='monotone'",
               column(3,sliderInput("LU_m_inc", "Bounds of Uniform Prior for increment between means of two adjacent DLs (m%)",
                                    min=0,max=200,step=5,value=c(0,50))),
             ),
           ),
           fluidRow(
             column(3,numericInput("n_pt_max","Max No. of Subjects per Cohort",6)),
             column(3,numericInput("n_pt_min","Min No. of Subjects to Evaluate FIX Posterior/Predicted Probability Criteria",3)),
             column(3,numericInput("dlt_thresh_r","Threshold of Posterior DLT Rate to Stop the Trial",0.15)),
             column(3,numericInput("dlt_thresh_c","Threshold of Posterior DLT Rate Probability to Stop the Trial",0.8)),
           ),
           fluidRow(
             column(3,numericInput("fix_thresh","Threshold of Predicted FIX Functional Activity > (Threshold of Predicted FIX) Probability to Stop the Trial",0.1)),
             column(3,sliderInput("fix_interval_ind", "Thresholds of Observed and Predicted (or Observed) FIX Functional Activity (%)",
                                  min=0,max=300,value=c(5,150),width="90%")),
             column(3,sliderInput("fix_interval_meanraw", "Thresholds of Posterior Mean of FIX Functional Activity (%)",
                                  min=0,max=300,value=c(35,70),width="90%")),
           )
          )
}

tab_conduct_study <- function() {
    tabPanel("Conduct Study",
             sidebarLayout(
                 sidebarPanel(
                     textInput("y_FIX", "Input FIX Functional Activity Data (y%) (Use ',' to Split Data)", 
                               "20, 60, 40, 10, 50, 25"),
                     textInput("x_DLT", "Input Binary DLT Data (Use ',' to Split Data)", 
                               "0,1,0,0,1,0"),
                     conditionalPanel(
                       condition="input.bayes_model=='same_cv' | input.bayes_model=='monotone'",
                       textInput("y_FIX_hist", "Input FIX Functional Activity Data (y%) before the current DL (Use ',' to Split Data)", 
                                 "2, 10, 5, 12, 9, 15, 20"),
                       textInput("DL_hist", "Input the corresponding DLs (Use ',' to Split Data)", 
                                 "1, 2, 2, 2, 3, 3, 3"),
                     ),
                     conditionalPanel(
                       condition="input.conduct_study=='Resulting Action'",
                       radioButtons("accelerated_titration", "Accelerated titration due to very low observed FIX of the 1st patient",
                                    c("Yes"=1,"No"=0),inline=TRUE),
                     )
                 ),
                 mainPanel(
                     tabsetPanel(type = "tabs", id = "conduct_study",
                                 tabPanel("Post Mean",
                                          h4("Kernel Smooth Density Plot of Posterior Mean of FIX Functional Activity"),
                                          plotOutput('density_m', width = "90%", height = "450px"),
                                          sliderInput("slider_m",
                                                      "Thresholds of Mean of FIX Functional Activity (%)",
                                                      min=0,max=200,step=5,value=c(35,70),width="90%"),
                                          sliderInput("slider_m_zoom", "Zoom in Mean of FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(0,300),width="90%"),
                                          ),
                                 tabPanel("Post Predict",
                                          h4("Kernel Smooth Density Plot of Predicted FIX Functional Activity"),
                                          plotOutput('density_y_tilde', width = "90%", height = "450px"),
                                          sliderInput("slider_y_tilde",
                                                      "Thresholds of Predicted FIX Functional Activity (%)",
                                                      min=0,max=200,step=5,value=c(5,150),width="90%"),
                                          sliderInput("slider_y_tilde_zoom",
                                                      "Zoom in Predicted FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(0,300),width="90%"),
                                          ),
                                 tabPanel("Post CV",
                                          h4("Kernel Smooth Density Plot of Posterior Coefficient of Variation of FIX Functional Activity"),
                                          plotOutput('density_cv', width = "90%", height = "450px"),
                                          sliderInput("slider_cv",
                                                      "Thresholds of Coefficient of Variation of FIX Functional Activity",
                                                      min=0,max=3,step=0.1,value=c(0.5,1.2),width="90%"),
                                          sliderInput("slider_cv_zoom",
                                                      "Zoom in Coefficient of Variation of FIX Functional Activity",
                                                      min=0,max=3,step=0.1,value=c(0,3),width="90%"),
                                          ),
                                 tabPanel("Summary Table of FIX",
                                          h4("Summary Table"),
                                          tableOutput("table_FIX"),
                                          ),
                                 tabPanel("Prior and Post DLT Rate",
                                          h4("Density Plot of Prior DLT Rate"),
                                          plotOutput('pdf_p_prior', width = "90%", height = "300px"),
                                          h4("Density Plot of Posterior DLT Rate"),
                                          plotOutput('pdf_p_post', width = "90%", height = "300px"),
                                          sliderInput("slider_p", "Threshold of DLT Rate",min=0,max=1,value=0.15,width="90%"),
                                          ),
                                 tabPanel("Resulting Action",
                                          h4("Action"),
                                          tableOutput('dose_action'),
                                 ),
                                 )
                     )
             )
             )
}

tab_simu <- function() {
    tabPanel("Simulation Results",
             sidebarLayout(
                 sidebarPanel(
                   radioButtons("bayes_model_select_SR", "Bayesian Model",
                                c("Independent DLs"="independent",
                                  "Same CV of FIX"="same_cv", 
                                  "Same CV and Monotone Increasing FIX"="monotone")),
                   radioButtons("m_true_SR", "True Mean of FIX Functional Activity (m%)",
                                c("20, 50, 80, 100" = 1, 
                                  "10, 20, 30, 80" = 2,
                                  "20, 60, 100, 120" = 3,
                                  "65, 85, 95, 120" = 4,
                                  "5, 20, 30, 45" = 5,
                                  "20, 35, 70, 100" = 6)),
                   radioButtons("cv_true_SR", "True Coefficient of Variation of FIX Functional Activity (CV)",
                                c("0.5" = 0.5, 
                                  "0.8" = 0.8, 
                                  "1.0" = 1, 
                                  "1.2" = 1.2)),
                   radioButtons("m_prior_SR", "Prior for Mean of FIX (m%)",
                                c("Uniform (20, 100)" = 1,
                                  "Uniform (0, 150)" = 2)),
                   radioButtons("cv_prior_SR", "Prior for Coefficient of Variation of FIX Functional Activity (CV)",
                                c("Uniform (0.4, 1.2)" = 1,
                                  "Uniform (0, 2)" = 2)),
                   strong("DLT Rate:"),
                   p("0.01, 0.01, 0.02, 0.02"),
                   strong("**Default values are used for the other design parameters."),
                 ),

                 mainPanel(
                   tabsetPanel(type = "tabs",
                               tabPanel("Scenario Plot",
                                        h4("Plot of DLT Rate and Mean of FIX Functional Activity"),
                                        plotOutput('plot_scenario_SR', width = "90%", height = "450px"),
                               ),
                               tabPanel("FIX Density",
                                        radioButtons("log_scale_SR", "FIX Functional Activity in Log-Scale",
                                                     c("No"= 0, "Yes"=1),inline=TRUE),
                                        conditionalPanel(
                                          condition="input.log_scale_SR==0",
                                          h4("Density Plot of FIX Functional Activity"),
                                          plotOutput('density_y_true_SR', width = "90%", height = "450px"),
                                          sliderInput("slider_y_true_SR", "Thresholds of FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(5,150),width="90%"),
                                          sliderInput("slider_y_true_SR_zoom", "Zoom in FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(0,300),width="90%"),
                                          sliderInput("slider_f_true_SR_zoom", "Zoom in Density Function",
                                                      min=0,max=0.1,step=0.01,value=0.05,width="90%"),
                                        ),
                                        conditionalPanel(
                                          condition="input.log_scale_SR==1",
                                          h4("Density Plot of Log FIX Functional Activity"),
                                          plotOutput('density_log_y_true_SR', width = "90%", height = "450px"), 
                                          sliderInput("slider_log_y_true_SR", "Thresholds of Log-Scale FIX Functional Activity",
                                                      min=0,max=7,step=0.1,value=c(1.5,5),width="90%"),
                                          sliderInput("slider_log_y_true_SR_zoom", "Zoom in Log-Scale FIX Functional Activity",
                                                      min=0,max=7,step=0.1,value=c(0,7),width="90%"),
                                          sliderInput("slider_f_log_true_SR_zoom", "Zoom in Density Function",
                                                      min=0,max=2,step=0.1,value=1,width="90%"),
                                        ),   
                               ),
                               tabPanel("Observed FIX Density",
                                        radioButtons("log_scale_observed_SR", "Observed FIX Functional Activity in Log-Scale",
                                                     c("No"= 0, "Yes"=1),inline=TRUE),
                                        conditionalPanel(
                                          condition="input.log_scale_observed_SR==0",
                                          h4("Density Plot of Observed FIX Functional Activity"),
                                          plotOutput('density_y_sim_SR', width = "90%", height = "450px"),
                                          sliderInput("slider_y_sim_SR", "Thresholds of FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(5,150),width="90%"),
                                          sliderInput("slider_y_sim_SR_zoom", "Zoom in FIX Functional Activity (%)",
                                                      min=0,max=300,step=5,value=c(0,300),width="90%"),
                                          sliderInput("slider_f_sim_SR_zoom", "Zoom in Density Function",
                                                      min=0,max=0.1,step=0.01,value=0.05,width="90%"),
                                        ),
                                        conditionalPanel(
                                          condition="input.log_scale_observed_SR==1",
                                          h4("Density Plot of Log Observed FIX Functional Activity"),
                                          plotOutput('density_log_y_sim_SR', width = "90%", height = "450px"),
                                          sliderInput("slider_log_y_sim_SR", "Thresholds of Log-Scale FIX Functional Activity",
                                                      min=0,max=7,step=0.1,value=c(1.5,5),width="90%"),
                                          sliderInput("slider_log_y_sim_SR_zoom", "Zoom in Log-Scale FIX Functional Activity",
                                                      min=0,max=7,step=0.1,value=c(0,7),width="90%"),
                                          sliderInput("slider_f_log_sim_SR_zoom", "Zoom in Density Function",
                                                      min=0,max=2,step=0.1,value=1,width="90%"),
                                        ),
                               ),
                               tabPanel("Operating Characteristics",
                                        h4("Summary Table of Operating Characteristics"),
                                        tableOutput('table_results_SR'),
                               ),
                   )
                )
             ),
            )
}

tab_simulation_study <- function(){
  tabPanel("Conduct Simulation",
           sidebarLayout(
             sidebarPanel(
               textInput("true_m", "True Mean of FIX Functional Activity for all Cohorts (y%) (Use ',' to Split)", "20, 50, 80, 100"),
               textInput("true_cv","True Coefficient of Variation (CV) of FIX Functional Activity for all Cohorts (Use ',' to Split for different CV)", 0.8),
               textInput("true_p", "True DLT Rate for all Cohorts (Use ',' to Split)", "0.01, 0.01, 0.02, 0.02"),
               textInput("ar_dose", "Upper and Lower Bounds of Ideal Admissible Doses (Use ',' to Split)", "2"),
               numericInput("n_sim","Number of Simulation Replicates",100),
               actionButton("button", "Run Simulation"),
               
             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Scenario Plot",
                                    h4("Plot of DLT Rate and Mean of FIX Functional Activity"),
                                    plotOutput('plot_scenario_CS', width = "90%", height = "450px"),
                           ),
                           tabPanel("FIX Density",
                                    radioButtons("log_scale_CS", "FIX Functional Activity in Log-Scale",
                                                 c("No"= 0, "Yes"=1),inline=TRUE),
                                    conditionalPanel(
                                      condition="input.log_scale_CS==0",
                                      h4("Density Plot of FIX Functional Activity"),
                                      plotOutput('density_y_true_CS', width = "90%", height = "450px"),
                                      sliderInput("slider_y_true_CS", "Thresholds of FIX Functional Activity (%)",
                                                  min=0,max=300,step=5,value=c(5,150),width="90%"),
                                      sliderInput("slider_y_true_CS_zoom", "Zoom in FIX Functional Activity",
                                                  min=0,max=300,step=5,value=c(0,300),width="90%"),
                                      sliderInput("slider_f_true_CS_zoom", "Zoom in Density Function",
                                                  min=0,max=0.1,step=0.01,value=0.05,width="90%"),
                                    ),
                                    conditionalPanel(
                                      condition="input.log_scale_CS==1",
                                      h4("Density Plot of Log FIX Functional Activity"),
                                      plotOutput('density_log_y_true_CS', width = "90%", height = "450px"),
                                      sliderInput("slider_log_y_true_CS", "Thresholds of Log-Scale FIX Functional Activity",
                                                  min=0,max=7,step=0.1,value=c(1.5,5),width="90%"),
                                      sliderInput("slider_log_y_true_CS_zoom", "Zoom in Log-Scale FIX Functional Activity",
                                                  min=0,max=7,step=0.1,value=c(0,7),width="90%"),
                                      sliderInput("slider_f_log_true_CS_zoom", "Zoom in Density Function",
                                                  min=0,max=2,step=0.1,value=1,width="90%"),
                                    ),   
                           ),
                           tabPanel("Observed FIX Density",
                                    radioButtons("log_scale_observed_CS", "Observed FIX Functional Activity in Log-Scale",
                                                 c("No"= 0, "Yes"=1),inline=TRUE),
                                    conditionalPanel(
                                      condition="input.log_scale_observed_CS==0",
                                      h4("Density Plot of Observed FIX Functional Activity"),
                                      plotOutput('density_y_sim_CS', width = "90%", height = "450px"),
                                      sliderInput("slider_y_sim_CS", "Thresholds of FIX Functional Activity (%)",
                                                  min=0,max=300,step=5,value=c(5,150),width="90%"),
                                      sliderInput("slider_y_sim_CS_zoom", "Zoom in FIX Functional Activity (%)",
                                                  min=0,max=300,step=5,value=c(0,300),width="90%"),
                                      sliderInput("slider_f_sim_CS_zoom", "Zoom in Density Function",
                                                  min=0,max=0.1,step=0.01,value=0.05,width="90%"),
                                    ),
                                    conditionalPanel(
                                      condition="input.log_scale_observed_CS==1",
                                      h4("Density Plot of Log Observed FIX Functional Activity"),
                                      plotOutput('density_log_y_sim_CS', width = "90%", height = "450px"),
                                      sliderInput("slider_log_y_sim_CS", "Thresholds of Log-Scale FIX Functional Activity",
                                                  min=0,max=7,step=0.1,value=c(1.5,5),width="90%"),
                                      sliderInput("slider_log_y_sim_CS_zoom", "Zoom in Log-Scale FIX Functional Activity",
                                                  min=0,max=7,step=0.1,value=c(0,7),width="90%"),
                                      sliderInput("slider_f_log_sim_CS_zoom", "Zoom in Density Function",
                                                  min=0,max=2,step=0.1,value=1,width="90%"),
                                    ),
                           ),
                           tabPanel("Operating Characteristics",
                                    h4("Summary Table of Operating Characteristics"),
                                    tableOutput('table_results_CS'),
                           ),
               )
             )
           )
  )
}

tab_other_settings <- function() {
  tabPanel("Other Settings",
           sidebarLayout(
             sidebarPanel(
               numericInput("seed","Choose a Seed for Simulation",100),
             ),
             
             mainPanel()
           )
  )
}


##-------------------------------------------------------------
##           DATA MANIPULATION: DEMO 1
##-------------------------------------------------------------


#-------------------------Conduct Study-----------------------

lst_design = reactive({
  list(n_pt_max = input$n_pt_max, 
       n_pt_min = input$n_pt_min, 
       dlt_thresh_r = input$dlt_thresh_r, 
       dlt_thresh_c = input$dlt_thresh_c,
       dlt_prior = c(input$alpha,input$beta), 
       fix_thresh = input$fix_thresh, 
       fix_prior_meanraw = input$LU_m, 
       fix_prior_cv = input$LU_cv, 
       fix_prior_meaninc = input$LU_m_inc,
       fix_interval_ind = input$fix_interval_ind, 
       fix_interval_meanraw = input$fix_interval_meanraw, 
       algorithm = desfix_algorithm_1,
       mean_raw = true_m(), 
       cv = true_cv(), 
       dlt_rates = true_p(), 
       ar_dose = ar_dose(),
       bayes_model = input$bayes_model,
       n_min_borrow = input$n_min_borrow)
})

# FIX
y_FIX = reactive({
    return(as.numeric(str_split(input$y_FIX,",", simplify = TRUE))) # percentage
})

data_hist = reactive({
  
  functional_activity = as.numeric(str_split(input$y_FIX_hist,",", simplify = TRUE))
  dose_level = as.numeric(str_split(input$DL_hist,",", simplify = TRUE))
  data_hist = data.frame(functional_activity,dose_level)
  
})


# Summary Table
table_FIX = reactive({
    table = data.frame(matrix(c(mean(y_FIX()),sd(y_FIX()),sd(y_FIX())/mean(y_FIX()),quantile(y_FIX())),nrow=1))
    colnames(table) = c("Mean (%)","SD (%)","CV","Min (%)","Q1 (%)","Median (%)","Q3 (%)","Max (%)")
    return(table)
})

# Posterior Sampling
fit_FIX = reactive({
    set.seed(input$seed)
    fit = desfix_bayes(y = y_FIX(),
                       L_m = input$LU_m[1],
                       U_m = input$LU_m[2],
                       L_cv = input$LU_cv[1],
                       U_cv = input$LU_cv[2],
                       L_m_inc = input$LU_m_inc[1],
                       U_m_inc = input$LU_m_inc[2],
                       data_hist = data_hist(),
                       n_min_borrow = input$n_min_borrow,
                       bayes_model = input$bayes_model)
    return(fit)
})

# Density plot
density_m = reactive({
  if (is.null(input$slider_m)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_m[,ncol(fit_FIX()$post_m)], # plot the current DL
                             cut = input$slider_m,
                             xlab_name = "Mean of FIX Functional Activity (m%)",
                             para_notation = bquote(m),percent = "Y",
                             from = 0,to = NA)
  return(pdf)
})

density_cv = reactive({
  if (is.null(input$slider_cv)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_cv[,ncol(fit_FIX()$post_cv)], # plot the current DL
                             cut = input$slider_cv,
                             xlab_name = "Coefficient of Variation of FIX Functional Activity (CV)",
                             para_notation = bquote(CV),percent = "N",
                             from = 0,to = NA)
  return(pdf)
})

density_y_tilde = reactive({
  if (is.null(input$slider_y_tilde)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_y_tilde[,ncol(fit_FIX()$post_y_tilde)], # plot the current DL
                             cut = input$slider_y_tilde,
                             xlab_name = bquote("Predicted FIX Functional Activity ("*tilde(y)*"%)"),
                             para_notation = bquote(tilde(y)),
                             percent = "Y",from = 0,to = NA)
  return(pdf)
})


# DLT
x_DLT = reactive({
  return(as.numeric(str_split(input$x_DLT,",", simplify = TRUE))) # percentage
})

# prior and poster density plot
pdf_p_prior = reactive({
  pdf = desfix_beta_binomial_desity(x = NA,
                                   n = NA,
                                   alpha = input$alpha,
                                   beta = input$beta,
                                   cut = input$slider_p,
                                   prior_post = "prior")
  return(pdf)
})

pdf_p_post = reactive({
  pdf = desfix_beta_binomial_desity(x = sum(x_DLT()),
                                   n = length(x_DLT()),
                                   alpha = input$alpha,
                                   beta = input$beta,
                                   cut = input$slider_p,
                                   prior_post = "post")
  return(pdf)
})

# Does Action
true_m = reactive({
  return(as.numeric(str_split(input$true_m,",", simplify = TRUE))) 
})

true_cv = reactive({
  return(as.numeric(str_split(input$true_cv,",", simplify = TRUE))) 
})

true_p = reactive({
  return(as.numeric(str_split(input$true_p,",", simplify = TRUE)))
})

ar_dose = reactive({
  return(as.numeric(str_split(input$ar_dose,",", simplify = TRUE)))
})


dose_action = reactive({
  
  #' Dose escalation algorithm
  #'
  #' 0:stay
  #'
  #' -12:stop based on predicted and posterior mean FIX given 6 subjects
  #' -13:stop based on predicted and posterior mean FIX given 3~5 subjects
  #' -14:stop due to both over posterior DLT given >=2 subjects and over observed FIX for any subject
  #' -15:stop due to over observed FIX for any subjects
  #' -16:stop due to over posterior DLT given >=2 subjects
  #' No stop rule for DLT when only 1 subject; Has stop rule for FIX for any subject.
  #'
  #' 11: escalate due to very low observed FIX of the 1st subjects
  #' 12: escalate based on predicted and posterior mean FIX given 6 subjects
  #' 13: escalate based on predicted and posterior mean FIX given 3~5 subjects
  #'
  #' 2:selected based on predicted and posterior mean FIX given 6 subjects
  
  set.seed(input$seed)
  
  cur_data = data.frame(y_FIX(),x_DLT())
  
  data_hist_with_AT = data_hist()
  if (input$accelerated_titration==1) {
    data_hist_with_AT$decision = 11 # Make sure accelerated titration. See function desfix_algorithm_1().
  } else if (input$accelerated_titration==0) {
    data_hist_with_AT$decision = 100 # Default setting. Posterior calculation is independent of decision. Set 100 in case being removed due to NA.
  }
  
  temp_dose_action = data.frame(Subj=1:length(y_FIX()),
                                Action=desfix_algorithm_1(cur_data = cur_data,
                                                          lst_design = lst_design(), 
                                                          data_hist = data_hist_with_AT))
  dose_action = data.frame(Subj=1:length(y_FIX()),Action=NA)
  dose_action$Action[temp_dose_action$Action==0] = paste("Stay")
  dose_action$Action[temp_dose_action$Action==-12] = paste("Stop based on predicted and posterior mean FIX given",input$n_pt_max,"subjects")
  dose_action$Action[temp_dose_action$Action==-13] = paste("Stop based on predicted and posterior mean FIX given",input$n_pt_min,"~",input$n_pt_max-1,"subjects")
  dose_action$Action[temp_dose_action$Action==-14] = paste("Stop due to both over posterior DLT given >= 2 subjects and over observed FIX for any subject")
  dose_action$Action[temp_dose_action$Action==-15] = paste("Stop due to over observed FIX for any subject")
  dose_action$Action[temp_dose_action$Action==-16] = paste("Stop due to over posterior DLT given >= 2 subjects")
  dose_action$Action[temp_dose_action$Action==11] = paste("Escalate due to very low observed FIX of the 1st patient")
  dose_action$Action[temp_dose_action$Action==12] = paste("Escalate based on predicted and posterior mean FIX given",input$n_pt_max,"subjects")
  dose_action$Action[temp_dose_action$Action==13] = paste("Escalate based on predicted and posterior mean FIX given",input$n_pt_min,"~",input$n_pt_max-1,"subjects")
  dose_action$Action[temp_dose_action$Action==2] = paste("Select")
  
  return(dose_action)
})

#---------------------Simulation Results------------------

# Values to map the data set
m_true = rbind(c(20, 50, 80, 100),
               c(10, 20, 30, 80),
               c(20, 60, 100, 120),
               c(65, 85, 95, 120),
               c(5, 20, 30, 45),
               c(20, 35, 70, 100))
m_prior = rbind(c(0,100),
                c(20,100),
                c(0,150),
                c(20,150))
cv_prior = rbind(c(0,2),
                 c(0.1,1.5),
                 c(0.4,1.2))

lst_design_SR = reactive({
  
  if (as.numeric(input$m_true_SR)==1) {
    ar_dose = 2
  } else if (as.numeric(input$m_true_SR)==2) {
    ar_dose = NA
  } else if (as.numeric(input$m_true_SR)==3) {
    ar_dose = 2
  } else if (as.numeric(input$m_true_SR)==4) {
    ar_dose = 1
  } else if (as.numeric(input$m_true_SR)==5) {
    ar_dose = 4
  } else if (as.numeric(input$m_true_SR)==6) {
    ar_dose = c(2,3)
  }
  list(n_pt_max = 6, 
       n_pt_min = 3, 
       dlt_thresh_r = 0.15, 
       dlt_thresh_c = 0.8,
       dlt_prior = c(0.1, 0.9), 
       fix_thresh = 0.1, 
       fix_prior_meanraw = m_prior[as.numeric(input$m_prior_SR),], 
       fix_prior_cv = cv_prior[as.numeric(input$cv_prior_SR),], 
       fix_prior_meaninc = NA,
       fix_interval_ind = c(5, 150), 
       fix_interval_meanraw = c(35, 70), 
       algorithm = desfix_algorithm_1,
       mean_raw = m_true[as.numeric(input$m_true_SR),], 
       cv = as.numeric(input$cv_true_SR), 
       dlt_rates = c(0.01,0.01,0.02,0.02), 
       ar_dose = ar_dose,
       bayes_model = input$bayes_model_select_SR,
       n_min_borrow = input$n_min_borrow)
  
})

# Saved simulation results

results_sim = reactive({
  results_sim = c()
  for (i in 1:50){
    results_sim = rbind(results_sim,get(load(paste0("Results/rst_",i,".Rdata"))))
  }
  return(results_sim)
})

results_sim_SR = reactive({

  
  # subset of data
  results_sim_SR = results_sim()[results_sim()$mean_raw    == input$m_true_SR & 
                                 results_sim()$cv          == input$cv_true_SR & 
                                 results_sim()$mean_pri    == input$m_prior_SR & 
                                 results_sim()$cv_pri      == input$cv_prior_SR &
                                 results_sim()$bayes_model == input$bayes_model_select_SR,]
  
  return(results_sim_SR)
})

# Plot of Scenario
plot_scenario_SR = reactive({
  return(desfix_plot_scenario(lst_design_SR()))
})

# Density Plot of true FIX for each DL
density_y_true_SR = reactive({
  return(lognormal_density_multiple(cut = input$slider_y_true_SR,
                                    true_m = m_true[as.numeric(input$m_true_SR),],
                                    true_cv = as.numeric(input$cv_true_SR),
                                    xlab_name = "FIX Functional Activity (y%)",
                                    para_notation = bquote(y), 
                                    logYN = "N", 
                                    percent = "Y"))
})
density_log_y_true_SR = reactive({
  return(lognormal_density_multiple(cut = input$slider_log_y_true_SR,
                                    true_m = m_true[as.numeric(input$m_true_SR),],
                                    true_cv = as.numeric(input$cv_true_SR),
                                    xlab_name = "FIX Functional Activity in Log-Scale (log(y))",
                                    para_notation = bquote(log(y)), 
                                    logYN = "Y", 
                                    percent = "N"))
})


# Density Plot of simulated FIX for each DL
density_y_sim_SR = reactive({
  pdf = desfix_kernel_density_multiple(cut = input$slider_y_sim_SR, 
                                       x = results_sim_SR()[,c("functional_activity","dose_level")],
                                       xlab_name = "FIX Functional Activity (y%)",
                                       para_notation = bquote(y),
                                       logYN = "N",
                                       percent = "Y",
                                       from = 0, to = NA)
  return(pdf)
})
density_log_y_sim_SR = reactive({
  pdf = desfix_kernel_density_multiple(cut = input$slider_log_y_sim_SR, 
                                       x = results_sim_SR()[,c("functional_activity","dose_level")],
                                       xlab_name = "FIX Functional Activity in Log-Scale (log(y))",
                                       para_notation = bquote(log(y)),
                                       logYN = "Y",
                                       percent = "N",
                                       from = NA, to = NA)
  return(pdf)
})

table_results_SR = reactive({
  return(desfix_summary(results_sim_SR(), lst_design_SR()))
})
  
#-------------------Conduct Simulation-----------------

results_sim_CS = eventReactive(input$button, {
  
  seed_v = input$seed
  lst_design_v = lst_design()
  cl = makeCluster(detectCores()-1)
  clusterEvalQ(cl, {library(dplyr); library(rstan)})
  clusterExport(cl, c("seed_v","lst_design_v","desfix_gen_data", "desfix_single_trial",
                      "desfix_generate_patients_cv", "desfix_generate_patients",
                      "get_prob_intervals","get_guarded","desfix_algorithm_1",
                      "desfix_single_trial","desfix_bayes"), 
                envir=environment())
  
  par_output_results = function(i) {
    set.seed(seed_v+i-1)
    dta_all = desfix_gen_data(lst_design_v)
    return(data.frame(desfix_single_trial(dta_all, lst_design_v),
                      rep=i,
                      cv=lst_design_v$cv,
                      mean_raw=rep(lst_design_v$mean_raw,each=lst_design_v$n_pt_max)))
  }
  temp_results_sim_CS = parSapply(cl, 1:input$n_sim, par_output_results)
  stopCluster(cl)
  
  results_sim_CS = c()
  for (i in 1:input$n_sim){
    results_sim_CS = rbind(results_sim_CS,data.frame(temp_results_sim_CS[,i]))
  }
  results_sim_CS = results_sim_CS[!is.na(results_sim_CS$decision),]
  return(results_sim_CS)
})


# Plot of Scenario for Simulation
plot_scenario_CS = reactive({
  return(desfix_plot_scenario(lst_design()))
})

# Density Plot of true FIX for each DL
density_y_true_CS = reactive({
  return(lognormal_density_multiple(cut = input$slider_y_true_CS,
                                    true_m = true_m(),
                                    true_cv = true_cv(),
                                    xlab_name = "FIX Functional Activity (y%)",
                                    para_notation = bquote(y), 
                                    logYN = "N", 
                                    percent = "Y"))
})
density_log_y_true_CS = reactive({
  return(lognormal_density_multiple(cut = input$slider_log_y_true_CS,
                                    true_m = true_m(),
                                    true_cv = true_cv(),
                                    xlab_name = "FIX Functional Activity in Log-Scale (log(y))",
                                    para_notation = bquote(log(y)), 
                                    logYN = "Y", 
                                    percent = "N"))
})


# Density Plot of simulated FIX for each DL
density_y_sim_CS = reactive({
  pdf = desfix_kernel_density_multiple(cut = input$slider_y_sim_CS, 
                                       x = results_sim_CS()[,c("functional_activity","dose_level")],
                                       xlab_name = "FIX Functional Activity (y%)",
                                       para_notation = bquote(y),
                                       logYN = "N",
                                       percent = "Y",
                                       from = 0, to = NA)
  return(pdf)
})
density_log_y_sim_CS = reactive({
  pdf = desfix_kernel_density_multiple(cut = input$slider_log_y_sim_CS,
                                       x = results_sim_CS()[,c("functional_activity","dose_level")],
                                       xlab_name = "FIX Functional Activity in Log-Scale (log(y))",
                                       para_notation = bquote(log(y)),
                                       logYN = "Y",
                                       percent = "N",
                                       from = NA, to = NA)
  return(pdf)
})


# Simulation Summary Table
table_results_CS = reactive({
  return(desfix_summary(results_sim_CS(), lst_design()))
})



  
  
