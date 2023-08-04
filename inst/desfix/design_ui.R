
##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------

##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_fix(),
                tab_dlt(),
                tab_simu(),
                tab_setting()
                )
}

tab_fix <- function() {
    tabPanel("FIX Functional Activity",
             sidebarLayout(
                 sidebarPanel(
                     textInput("y", "Input FIX Functional Activity Data (y%) (Use ',' to Split Data)", "20, 60, 40, 10, 50, 25"),
                     sliderInput("LU_m", "Bounds of Uniform Prior for Mean (m%)",min=0,max=500,step=5,value=c(0,300)),
                     sliderInput("LU_cv", "Bounds of Uniform Prior for Coefficient of Variation (cv)",min=0,max=5,step=0.1,value=c(0.4,2)),
                     ),
                 mainPanel(
                     tabsetPanel(type = "tabs",
                                 tabPanel("Post Mean",
                                          h4("Kernel Smooth Density Plot of Posterior Mean of FIX Functional Activity"),
                                          plotOutput('density_m', width = "90%", height = "450px"),
                                        # h4("Histgram of Posterior Sample of Mean of FIX"),
                                        # plotOutput('Hist_m', width = "90%", height = "450px"),
                                          htmlOutput("slider_m"),
                                          htmlOutput("slider_m_zoom"),
                                          ),
                                 tabPanel("Post Predict",
                                          h4("Kernel Smooth Density Plot of Predicted FIX Functional Activity"),
                                          plotOutput('density_y_tilde', width = "90%", height = "450px"),
                                        # h4("Histgram of Posterior Sample of Predicted FIX Functional Activity"),
                                        # plotOutput('Hist_y_tilde', width = "90%", height = "450px"),
                                          htmlOutput("slider_y_tilde"),
                                          htmlOutput("slider_y_tilde_zoom"),
                                          ),
                                 tabPanel("Post CV",
                                          h4("Kernel Smooth Density Plot of Posterior Coefficient of Variation of FIX Functional Activity"),
                                          plotOutput('density_cv', width = "90%", height = "450px"),
                                        # h4("Histgram of Posterior Sample of Coefficient of Variation of FIX Functional Activity"),
                                        # plotOutput('Hist_cv', width = "90%", height = "450px"),
                                          htmlOutput("slider_cv"),
                                          htmlOutput("slider_cv_zoom"),
                                          ),
                                 tabPanel("Summary Table",
                                          h4("Summary Table"),
                                          tableOutput("table_FIX"),
                                          )
                                 )
                     )
             )
             )
}

tab_dlt <- function() {
    tabPanel("DLT Rate",
             sidebarLayout(
                 sidebarPanel(
                     numericInput("x", "Input No. of Responses", 2),
                     numericInput("n", "Input No. of Subjects", 6),
                     numericInput("alpha","Parameter alpha of Beta Prior for DLT Rate",0.1),
                     numericInput("beta","Parameter beta of Beta Prior for DLT Rate",0.9),
                     ),

                 mainPanel(
                     h4("Density Plot of Prior DLT Rate"),
                     plotOutput('pdf_p_prior', width = "90%", height = "300px"),
                     h4("Density Plot of Posterior DLT Rate"),
                     plotOutput('pdf_p_post', width = "90%", height = "300px"),
                     sliderInput("slider_p", "Threshold of DLT Rate",min=0,max=1,value=0.15,width="90%"),
                     )
             ))

}

tab_setting <- function() {
    tabPanel("Other Settings",
             sidebarLayout(
                 sidebarPanel(
                     numericInput("seed","Choose a Seed for Simulation",100),
                     ),

                 mainPanel()
             )
             )
}

tab_simu <- function() {
    tabPanel("Operating Characteristics",
             sidebarLayout(
                 sidebarPanel(
                     radioButtons("cv_true_sc", "True Coefficient of Variation of FIX Functional Activity (cv)",
                                  c("0.5" = 0.5, "0.8" = 0.8, "1.0" = 1, "1.2" = 1.2)),
                     radioButtons("m_prior_sc", "Prior for Mean of FIX (m)",
                                  c("Uniform (0%, 300%)" = 1, "Uniform (10%, 100%)" = 2)),
                     radioButtons("cv_prior_sc", "Prior for Coefficient of Variation of FIX Functional Activity (cv)",
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


##-------------------------------------------------------------
##           DATA MANIPULATION: DEMO 1
##-------------------------------------------------------------



results = c()
for (i in 1:50){
    results = rbind(results,get(load(paste0("Results/rst_",i,".Rdata"))))
}

DL = unique(results$dose_level)
DLT_rate = c(1,2,3,5)
m_true = c(20,50,80,100)
m_prior = rbind(c(0,300),c(10,100))
cv_prior = rbind(c(0,2),c(0.1,1.5))

                                        #------------------------FIX--------------------------
y = reactive({
    return(as.numeric(str_split(input$y,",", simplify = TRUE))) # percentage
})

                                        #-----------Summary Table-----------
table_FIX = reactive({
    table = data.frame(matrix(c(mean(y()),sd(y()),sd(y())/mean(y()),quantile(y())),nrow=1))
    colnames(table) = c("Mean (%)","SD (%)","CV","Min (%)","Q1 (%)","Median (%)","Q3 (%)","Max (%)")
    return(table)
})

                                        #----------Posterior Sampling----------
fit_FIX = reactive({
    set.seed(input$seed)
    fit = FIX_Bayesian(y(),input$LU_m[1],input$LU_m[2],input$LU_cv[1],input$LU_cv[2])
    return(fit)
})









  density_m = reactive({
    if (is.null(input$slider_m)){return(0)} # temperal value before having input value
    pdf = fun_kernel_density(x = fit_FIX()$post_m,
                             cut = input$slider_m,
                             xlab_name = "Mean of FIX Functional Activity (m%)",
                             para_notation = bquote(m),percent = "Y",
                             from = 0,to = NA)
    return(pdf)
  })

    density_cv = reactive({
    if (is.null(input$slider_cv)){return(0)} # temperal value before having input value
    pdf = fun_kernel_density(x = fit_FIX()$post_cv,
                             cut = input$slider_cv,
                             xlab_name = "Coefficient of Variation of FIX Functional Activity (cv)",
                             para_notation = bquote(cv),percent = "N",
                             from = 0,to = NA)
    return(pdf)
  })

    density_y_tilde = reactive({
    if (is.null(input$slider_y_tilde)){return(0)} # temperal value before having input value
    pdf = fun_kernel_density(x = fit_FIX()$post_y_tilde,
                             cut = input$slider_y_tilde,
                             xlab_name = bquote("Predicted FIX Functional Activity ("*tilde(y)*"%)"),
                             para_notation = bquote(tilde(y)),
                             percent = "Y",from = 0,to = NA)
    return(pdf)
  })


  #------------------------DLT Rate--------------------------

  pdf_p_prior = reactive({
    pdf = fun_beta_binomial_desity(x = NA,
                                   n = NA,
                                   alpha = input$alpha,
                                   beta = input$beta,
                                   cut = input$slider_p,
                                   prior_post = "prior")
    return(pdf)
  })

  pdf_p_post = reactive({
    pdf = fun_beta_binomial_desity(x = input$x,
                                   n = input$n,
                                   alpha = input$alpha,
                                   beta = input$beta,
                                   cut = input$slider_p,
                                   prior_post = "post")
    return(pdf)
  })

  #----------------Operating Characteristics------------

  plot_scenario = reactive({

  })
