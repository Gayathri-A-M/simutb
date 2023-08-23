
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
           h4("Hyper Parameters for Prior Distributions"),
           fluidRow(
             column(3,sliderInput("LU_m", "Bounds of Uniform Prior for Mean (m%)",min=0,max=500,step=5,value=c(0,300))),
             column(3,sliderInput("LU_cv", "Bounds of Uniform Prior for Coefficient of Variation (CV)",min=0,max=5,step=0.1,value=c(0.4,2))),
             column(3,numericInput("alpha","Parameter alpha of Beta Prior for DLT Rate",0.1)),
             column(3,numericInput("beta","Parameter beta of Beta Prior for DLT Rate",0.9)),
           ),

           h4("Parameters for Study Design"),
           fluidRow(
             column(3,numericInput("n_pt_max","Max No. of Subjects per Cohort",6)),
             column(3,numericInput("n_pt_min","Min No. of Subjects to Evaluate Posterior Probability Criteria",3)),
             column(3,numericInput("dlt_thresh_r","Threshold of DLT Rate to Stop the Trial",0.15)),
             column(3,numericInput("dlt_thresh_c","Threshold of DLT Rate Probability to Stop the Trial (%)",80)),
           ),
           fluidRow(
             column(3,numericInput("fix_thresh","Threshold of FIX Functional Activity Probability to Stop the Trial (%)",10)),
             column(3,sliderInput("fix_interval_ind", "Threshold of Predicted FIX Functional Activity to Continue the Trial (%)",
                                  min=0,max=300,value=c(5,150),width="90%")),
             column(3,sliderInput("fix_interval_meanraw", "Threshold of Mean of FIX Functional Activity to Continue the Trial",
                                  min=0,max=300,value=c(35,70),width="90%")),
           )
          )
}

tab_conduct_study <- function() {
    tabPanel("Conduct Study",
             sidebarLayout(
                 sidebarPanel(
                     textInput("y_FIX", "Input FIX Functional Activity Data (y%) (Use ',' to Split Data)", "20, 60, 40, 10, 50, 25"),
                     textInput("x_DLT", "Input Binary DLT Data (Use ',' to Split Data)", "0,1,0,0,1,0"),
                     ),
                 mainPanel(
                     tabsetPanel(type = "tabs",
                                 tabPanel("Post Mean",
                                          h4("Kernel Smooth Density Plot of Posterior Mean of FIX Functional Activity"),
                                          plotOutput('density_m', width = "90%", height = "450px"),
                                          htmlOutput("slider_m"),
                                          htmlOutput("slider_m_zoom"),
                                          ),
                                 tabPanel("Post Predict",
                                          h4("Kernel Smooth Density Plot of Predicted FIX Functional Activity"),
                                          plotOutput('density_y_tilde', width = "90%", height = "450px"),
                                          htmlOutput("slider_y_tilde"),
                                          htmlOutput("slider_y_tilde_zoom"),
                                          ),
                                 tabPanel("Post CV",
                                          h4("Kernel Smooth Density Plot of Posterior Coefficient of Variation of FIX Functional Activity"),
                                          plotOutput('density_cv', width = "90%", height = "450px"),
                                          htmlOutput("slider_cv"),
                                          htmlOutput("slider_cv_zoom"),
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
    tabPanel("Operating Characteristics",
             sidebarLayout(
                 sidebarPanel(
                     radioButtons("cv_true_sc", "True Coefficient of Variation of FIX Functional Activity (CV)",
                                  c("0.5" = 0.5, "0.8" = 0.8, "1.0" = 1, "1.2" = 1.2)),
                     radioButtons("m_prior_sc", "Prior for Mean of FIX (m)",
                                  c("Uniform (0%, 300%)" = 1, "Uniform (10%, 100%)" = 2)),
                     radioButtons("cv_prior_sc", "Prior for Coefficient of Variation of FIX Functional Activity (CV)",
                                  c("Uniform (0, 2)" = 1, "Uniform (0.1, 1.5)" = 2)),
                     ),

                 mainPanel(
                     h4("Plot of DLT Rate and Mean of FIX Functional Activity"),
                     plotOutput('plot_scenario', width = "90%", height = "450px"),
                     h4("Summary Table of Operating Characteristics"),
                     tableOutput('table_results'),
                     )
             ),
             )
}

tab_simulation_study <- function(){
  tabPanel("Simulation Study",
           sidebarLayout(
             sidebarPanel(
               textInput("true_m", "True Mean of FIX Functional Activity for all Cohorts (y%) (Use ',' to Split Data)", "20, 50, 80, 100"),
               numericInput("true_cv","True Coefficient of Variation (CV) of FIX Functional Activity for all Cohorts", 0.5),
               textInput("true_p", "True DLT Rate for all Cohorts (x%) (Use ',' to Split Data)", "1, 2, 3, 5"),
               textInput("ar_dose", "Upper and Lower Bounds of Ideal Admissible Doses (Use ',' to Split)", "2"),
               numericInput("n_sim","Number of Simulation Replicates",100),
               actionButton("button", "Run Simulation"),

             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Scenario Plot",
                                    h4("Plot of DLT Rate and Mean of FIX Functional Activity"),
                                    plotOutput('plot_scenario_sim', width = "90%", height = "450px"),
                           ),
                           tabPanel("Observed FIX Density",
                                    h4("Density Plot of Observed FIX Functional Activity"),
                                    plotOutput('density_y_sim', width = "90%", height = "450px"),
                           ),
                           tabPanel("Observed Log FIX Density",
                                    h4("Density Plot of Log Observed FIX Functional Activity"),
                                    plotOutput('density_log_y_sim', width = "90%", height = "450px"),
                           ),
                           tabPanel("Operating Characteristics",
                                    h4("Summary Table of Operating Characteristics"),
                                    tableOutput('table_results_sim'),
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

results = c()
for (i in 1:50){
    results = rbind(results,get(load(paste0("Results/rst_",i,".Rdata"))))
}

DL = unique(results$dose_level)
DLT_rate   = c(1,2,3,5)
m_true     = c(20,50,80,100)
m_prior    = rbind(c(0,300),c(10,100))
cv_prior   = rbind(c(0,2),c(0.1,1.5))

#-------------------------Conduct Study-----------------------

# FIX
y_FIX = reactive({
    return(as.numeric(str_split(input$y_FIX,",", simplify = TRUE))) # percentage
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
    fit = desfix_bayes(y_FIX(),input$LU_m[1],input$LU_m[2],input$LU_cv[1],input$LU_cv[2])
    return(fit)
})

# Density plot
density_m = reactive({
  if (is.null(input$slider_m)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_m,
                             cut = input$slider_m,
                             xlab_name = "Mean of FIX Functional Activity (m%)",
                             para_notation = bquote(m),percent = "Y",
                             from = 0,to = NA)
  return(pdf)
})

density_cv = reactive({
  if (is.null(input$slider_cv)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_cv,
                             cut = input$slider_cv,
                             xlab_name = "Coefficient of Variation of FIX Functional Activity (CV)",
                             para_notation = bquote(CV),percent = "N",
                             from = 0,to = NA)
  return(pdf)
})

density_y_tilde = reactive({
  if (is.null(input$slider_y_tilde)){return(0)} # temperal value before having input value
  pdf = desfix_kernel_density(x = fit_FIX()$post_y_tilde,
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
lst_design = reactive({
    list(n_pt_max = input$n_pt_max,
         n_pt_min = input$n_pt_min,
         dlt_thresh_r = input$dlt_thresh_r,
         dlt_thresh_c = input$dlt_thresh_c,
         dlt_prior = c(input$alpha,input$beta),
         fix_thresh = input$fix_thresh,
         fix_prior_meanraw = input$LU_m,
         fix_prior_cv = input$LU_cv,
         fix_interval_ind = input$fix_interval_ind,
         fix_interval_meanraw = input$fix_interval_meanraw,
         algorithm = desfix_algorithm_1,
         mean_raw = true_m(),
         cv = input$true_cv,
         dlt_rates = true_p(),
         ar_dose = ar_dose())
})

dose_action = reactive({

  #' -12:stop after 6 patients
  #' -13:stop before 6 patients
  #' -14:stop for both DLT and Over FIX
  #' -15:stop for Over FIX
  #' -16:stop for DLT
  #'
  #' 11: escalate after the 1st patient
  #' 12: escalate after 6 patients
  #' 13: escalate before 6 patients
  #'
  #' 2:selected

  cur_data = data.frame(y_FIX(),x_DLT())
  temp_dose_action = data.frame(Subj=1:length(y_FIX()),Action=desfix_algorithm_1(cur_data, lst_design()))
  dose_action = data.frame(Subj=1:length(y_FIX()),Action=NA)
  dose_action$Action[temp_dose_action$Action==0] = "Stay"
  dose_action$Action[temp_dose_action$Action==-12] = "Stop after 6 patients"
  dose_action$Action[temp_dose_action$Action==-13] = "Stop before 6 patients"
  dose_action$Action[temp_dose_action$Action==-14] = "Stop for both DLT and over FIX"
  dose_action$Action[temp_dose_action$Action==-15] = "Stop for over FIX only"
  dose_action$Action[temp_dose_action$Action==-16] = "Stop for DLT only"
  dose_action$Action[temp_dose_action$Action==11] = "Escalate after the 1st patient"
  dose_action$Action[temp_dose_action$Action==12] = "Escalate after 6 patients"
  dose_action$Action[temp_dose_action$Action==13] = "Escalate before 6 patients"
  dose_action$Action[temp_dose_action$Action==2] = "Select"

  return(dose_action)
})

#---------------------Operating Characteristics------------------

plot_scenario = reactive({
  AR_location_L = 2
  AR_location_U = 2

  data_plot = data.frame(Prob=c(DLT_rate,m_true,rep(35,length(DLT_rate)),rep(70,length(DLT_rate))),
                         DL=DL,Type=factor(rep(1:4,each=length(DL))))
  levels(data_plot$Type) = c("DLT Rate (%)","Mean of FIX F. A. (%)",
                             "Lower Bound of Target Mean of FIX F. A. (%)",
                               "Upper Bound of Target Mean of FIX F. A. (%)")

  plot_scenario = ggplot(data_plot,aes(x=DL,y=Prob,group=Type)) +
    geom_line(aes(color=Type,linetype=Type),size=1) +
    geom_point(aes(shape=Type,color=Type),size=4) +
    geom_rect(aes(xmin=AR_location_L-0.2,
                  xmax=AR_location_U+0.2,
                  ymin=0,ymax=100,
                  fill=""),
              alpha=0,color="green3") +
    scale_shape_manual(values=c(19, 19, NA, NA)) +
    scale_color_manual(values=c("red","deepskyblue2","grey","black")) +
    scale_fill_manual(values = "green3",labels="Ideal Dose") +
    ylim(0,100) + ylab("Probability") + xlab("Dose Level") +
    theme_bw() + theme(text=element_text(size=20),legend.position="top",
                       legend.title=element_blank(),
                       legend.key.height= unit(1, 'cm'),
                       legend.key.width= unit(2, 'cm')) +
    guides(fill=guide_legend(ncol=1),color=guide_legend(ncol=1),
           shape=guide_legend(ncol=1))

  return(plot_scenario)
})

table_results = reactive({

  results$dose_level = factor(results$dose_level,levels=1:length(m_true))

  table_results = c()
  sub_results = results[results$cv==input$cv_true_sc & results$mean_pri==input$m_prior_sc & results$cv_pri==input$cv_prior_sc,]
  n_rep = max(sub_results$rep)
  prob_sel = table(sub_results[sub_results$decision==2,c("dose_level")])/n_rep
  prob_stop = sum(sub_results$decision==-1)/n_rep
  prob_more_h = sum(sub_results$decision==1 & sub_results$dose_level==4)/n_rep
  ave_n = apply(table(sub_results[,c("dose_level","rep")]),1,mean)
  ave_n_total = mean(apply(table(sub_results[,c("dose_level","rep")]),2,sum))
  table_results = rbind(table_results,
                        data.frame(Cohort=1:length(m_true),
                                   DLT_rate=DLT_rate,
                                   m_true=m_true,
                                   prob_sel=as.numeric(round(prob_sel*100,1)),
                                   ave_n=round(ave_n,2),
                                   ave_n_total=c(rep("",length(m_true)-1),round(ave_n_total,2)),
                                   prob_stop = c(rep("",length(m_true)-1),round(prob_stop*100,1)),
                                   prob_more_h=c(rep("",length(m_true)-1),round(prob_more_h*100,1))))
  colnames(table_results) = c("Cohort","DLT Rate (%)","Mean of FIX F. A. (%)",
                              "Prob. Sel. (%)","Ave. No. of Subj. per Cohort",
                              "Ave. Total No. of Subj.","Prob. Stop (%)","Prob. Sel.>4 (%)")
  return(table_results)
})

#-------------------Simulation Study-----------------

true_m = reactive({
  return(as.numeric(str_split(input$true_m,",", simplify = TRUE)))
})

true_p = reactive({
  return(as.numeric(str_split(input$true_p,",", simplify = TRUE)))
})

ar_dose = reactive({
  return(as.numeric(str_split(input$ar_dose,",", simplify = TRUE)))
})

results_sim = eventReactive(input$button, {

  ## seed_v       = input$seed
    ## cl = makeCluster(detectCores()-1)
  ## clusterExport(cl, c("seed_v","lst_design_v","desfix_gen_data", "desfix_single_trial",
  ##                     "desfix_generate_patients_cv", "desfix_generate_patients",
  ##                     "get_prob_intervals","get_guarded","desfix_algorithm_1",
  ##                     "desfix_single_trial","desfix_bayes","stan","extract"),
  ##               envir=environment())

  ## par_output_results = function(i) {
  ##   set.seed(seed_v+i-1)
  ##   dta_all = desfix_gen_data(lst_design_v)
  ##   return(data.frame(desfix_single_trial(dta_all, lst_design_v),rep=i,cv=lst_design_v$cv,
  ##                                         mean_raw=rep(lst_design_v$mean_raw,each=lst_design_v$n_pt_max)))
  ## }
  ## temp_results_sim = parSapply(cl, 1:input$n_sim, par_output_results)
  ## stopCluster(cl)

  ## results_sim = c()
  ## for (i in 1:input$n_sim){
  ##   results_sim = rbind(results_sim,data.frame(temp_results_sim[,i]))
  ## }
  ## results_sim = results_sim[!is.na(results_sim$decision),]
  ## return(results_sim)

  ## simulation study
  xx           <- stb_create_design("dose_fix")
  stb_para(xx) <- lst_design()
  zz           <- stb_create_simustudy(xx,
                                       n_rep    = input$n_sim,
                                       n_core   = 5,
                                       seed     = input$seed,
                                       save_raw = TRUE,
                                       refresh  = 0,
                                       cores    = 1)

  zz
})


# Plot of Scenario for Simulation
plot_scenario_sim = reactive({
    xx           <- stb_create_design("dose_fix")
    stb_para(xx) <- lst_design()
    stb_plot_design(xx)
})


# Density Plot of simulated FIX for each DL
density_y_sim = eventReactive(input$button, {
  pdf = desfix_kernel_density_multiple(x = results_sim()[,c("functional_activity","dose_level")],
                                       xlab_name = "Mean of FIX Functional Activity (m%)",
                                       para_notation = bquote(m),percent = "Y",
                                       from = 0, to = NA, logYN = "N")
  return(pdf)
})

density_log_y_sim = eventReactive(input$button, {
  pdf = desfix_kernel_density_multiple(x = results_sim()[,c("functional_activity","dose_level")],
                                       xlab_name = "Mean of FIX Functional Activity (m%)",
                                       para_notation = bquote(m),percent = "Y",
                                       from = NA, to = NA, logYN = "Y")
  return(pdf)
})


# Simulation Summary Table
table_results_sim = eventReactive(input$button, {
  return(desfix_summary(results_sim(), lst_design()))
})
