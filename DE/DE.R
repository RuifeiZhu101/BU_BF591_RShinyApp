## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)
library(tidyverse)
source("de_functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Differential Expression"),
    p("Differential expression identifies which genes, if any, 
      are implicated in a specific biological comparison. 
      This component allows the user to load and explore a differential expression dataset."),
    sidebarLayout(
        sidebarPanel(
            fileInput("de_file",
                      "Load a differential expression dataset in CSV format:",
                      placeholder = "DEseq2_diffexp.csv",
                      accept = ".csv"),
            p("A volcano plot can be generated with", strong('"log2 fold-change"'),"on the x-axis and",
              strong('"p-adjusted"'),"on the y-axis."),
            
            radioButtons(inputId = 'xaxis', label = 'Choose the column for the x-axis', choices = c('baseMean',
                                                                                                    'log2FoldChange', 
                                                                                                    'lfcSE',
                                                                                                    'stat',
                                                                                                    'pvalue',
                                                                                                    'padj'), 
                         selected = 'log2FoldChange'),
            radioButtons(inputId = 'yaxis', label = 'Choose the column for the x-axis', choices = c('baseMean',
                                                                                                    'log2FoldChange', 
                                                                                                    'lfcSE',
                                                                                                    'stat',
                                                                                                    'pvalue',
                                                                                                    'padj'), 
                         selected = 'padj'),
            
            # Input 4: Color input
            colourInput("base",
                        "Base point color",
                        "#22577A"),
            
            # Input 5: Color input
            colourInput("highlight",
                        "Highlight point color",
                        "#FFCF56"),
            
            # Input 6: Slider bar to select the magnitude of adjusted p-value to color on the graph.
            sliderInput("slider", 
                        p("Select the magnitude of the p adjusted coloring:"),
                        min = -30, max = 0, value = -10),
            # Submit botton
            actionButton("submit_de",
                         "Plot", 
                         icon = icon("redo"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Plot",
                     plotOutput("volcano")
            ),
            tabPanel("Table",
                     dataTableOutput("de_table"))
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #' load_data on submission
  load_de_data <- reactive({
    req(input$de_file)
    deres_df <- read.csv(input$de_file$datapath) 
  return(deres_df)})
  
  #render volcano
  output$volcano <- renderPlot({
    req(input$submit_de)
    deseq_res <- load_de_data()
    volcano_plot(deseq_res,
                 input$xaxis, 
                 input$yaxis, 
                 input$slider,
                 input$base, 
                 input$highlight)
    })
  
  # render data table
  output$de_table <- renderDataTable({
    draw_table(deseq_res,input$slider)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
