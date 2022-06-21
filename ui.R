# User interface for simulation of non-differential misclassification of disease under (near)-perfect specificity and imperfect yet (expected) non-differential sensitivity

library(shiny)

shinyUI(fluidPage(

    # Application title
    titlePanel("NDMD Simulation"),

    # Sidebar with inputs for simulation parameters
    sidebarLayout(
        sidebarPanel(
          h3("Simulation parameters"),
          hr(),
          radioButtons("niter", 
                       "Number of iterations:", 
                       choices=list("100"=1e2, "1,000"=1e3, "10,000"=1e4, "100,000"=1e5, "1,000,000"=1e6),
                       selected=1e4),
          sliderInput("ip_set",
                      "Risk in the unexposed:",
                      min = 0.01,
                      max = 0.50,
                      value = 0.05),
          sliderInput("pe_set",
                      "Exposure prevalence:",
                      min=0.05,
                      max=0.50,
                      value=0.10),
          numericInput("r_set", 
                       "True risk ratio:",
                       value=1),
          sliderInput("se_set",
                      "Sensitivity of disease classification:",
                      min=0.5,
                      max=1.0,
                      value=0.8),
          sliderInput("sp_set",
                      "Specificity of disease classification:",
                      min=0.90,
                      max=1.0,
                      value=1.0),
          actionButton("run", "Run simulation")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            #plotOutput("distPlot")
        )
    )
))
