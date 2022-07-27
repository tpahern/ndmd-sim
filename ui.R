# User interface for simulation of non-differential misclassification of disease under (near)-perfect specificity and imperfect yet (expected) non-differential sensitivity

library(shiny)
library(shinythemes)
library(shinybusy)

shinyUI(fluidPage(theme=shinytheme("spacelab"),

    # Application title
    titlePanel("Simulation of non-differential misclassification of disease under (near-) perfect specificity",
               windowTitle = "NDMD Simulation"),
    br(),
    
    # Busy indicator
    # add_busy_spinner(spin="scaling-squares",
    #                 position="full-page"
    #                 ),
    
    # Smails busy GIF
    # add_busy_gif(src = "https://media.giphy.com/media/HQRgg6ks7nkyY/giphy.gif", height = 70, width = 70, position="full-page"),
    
    # Call me Al busy GIF
    # add_busy_gif(src = "https://media.giphy.com/media/10oguaVKJJqPZe/giphy.gif", height = 100, width = 100, position="full-page"),
    
    # A-Team busy GIF
    add_busy_gif(src = "https://media.giphy.com/media/pvwWuq1gdckNy/giphy.gif", height = 100, width = 100, position="full-page"),
    
    
    # Sidebar with inputs for simulation parameters
    sidebarLayout(
        sidebarPanel(width=3,
          h3("Inputs"),
          hr(),
          radioButtons("randseed",
                       "Seed type",
                       choices=list("Random (system clock)"=1,
                                    "Fixed (reproducible simulations)"=0)),
          radioButtons("niter", 
                       "Number of iterations", 
                       choices=list("10,000"=1e4, 
                                    "100,000"=1e5, 
                                    "1,000,000"=1e6)),
          radioButtons("studyn",
                      "Study size",
                      choices=list("10,000"=1e4,
                                   "50,000"=5e4,
                                   "100,000"=1e5,
                                   "500,000"=5e5,
                                   "1,000,000"=1e6)),
          sliderInput("ip_set",
                      "Risk in the unexposed",
                      min = 0.01,
                      max = 0.50,
                      value = 0.05),
          sliderInput("pe_set",
                      "Exposure prevalence",
                      min=0.05,
                      max=0.50,
                      step=0.05,
                      value=0.10),
          sliderInput("r_set", 
                      "True risk ratio",
                      min=1.0,
                      max=5.0,
                      step=0.1,
                      value=1.5
                      ),
          sliderInput("se_set",
                      "Sensitivity of disease classification",
                      min=0.5,
                      max=1.0,
                      value=0.8),
          sliderInput("sp_set",
                      "Specificity of disease classification",
                      min=0.90,
                      max=1.0,
                      value=0.99),
          hr(),
        
          actionButton("run", "Run simulation", class="btn-block",
                       icon = icon("thumbs-up", 
                                   lib="font-awesome")),
          br(),
          #br(),
        
          conditionalPanel(condition = "input.run > 0",
            downloadButton("downloadData", "Download simulation dataset", class="btn-block",
                         icon=icon("file-csv", 
                                   lib="font-awesome"))),
          br(),
          #br(),
          
          a(shiny::actionButton(inputId = "email1", 
                                label = "Contact author", 
                                class = "btn-block",
                                icon = icon("envelope", lib = "font-awesome")),
            href="mailto:02tahern@med.uvm.edu")
        ),

        # Output display
        mainPanel(

          conditionalPanel(
          condition = "input.run > 0",
          h3("Bias parameter distributions"),
          fluidRow(
            column(6,
                   h4("Sensitivity"),
                   plotOutput("se_jointdist"),
                   ),
            column(6,
                   h4("Specificity"),
                   plotOutput("sp_jointdist"),
                   ),
            ),

          hr(),
          
          h3("Simulation results"),
          fluidRow(
            column(6,
                   h4("Distribution of simulated true risk ratios"),
                   plotOutput("dens_iprc"),
                   verbatimTextOutput("iprc_sum"),
                   ),
            column(6,
                   h4("Distribution of misclassified risk ratios"),
                   plotOutput("dens_iprms"),
                   verbatimTextOutput("iprm_s_sum"),
                   verbatimTextOutput("iprm_areas")
                   )
            ),
          fluidRow(
            column(12,
                   h4("Bias factor as a function of departure from non-differential sensitivity and specificity (\u0394)"),
                   plotOutput("deltascat", height="400px")
                   )
            )
          )
        )
    )
    )
    )
