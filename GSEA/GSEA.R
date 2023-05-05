## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)
library(tidyverse)
source("GSEA_functions.R")


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Gene Set Enrichment Analysis with fgsea"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fgsea_file",
                "Upload your own fgsea results data in CSV format:",
                placeholder = "fgsea_res.csv",
                accept = ".csv"),
      actionButton("submit_fgsea","Upload"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Top Pathways",
                 sidebarLayout(
                   sidebarPanel(sliderInput("top_n", "Number of top pathways to plot by adjusted p-value:",
                                             min = 1, max = 50, value = 10),
                                 tags$hr(),
                                 tags$p("Click on a bar to display table entry"),
                                 verbatimTextOutput("selected_pathway")),
                   mainPanel(
                     plotOutput("pathways_barplot")
                   )
                 )),
        tabPanel("GSEA Table",
                 sidebarLayout(
                   sidebarPanel(
                     sliderInput("padj_cutoff", "Filter table by adjusted p-value (padj):",
                                 min = 0, max = 0.1, step = 0.01, value = 0.05),
                     radioButtons("nes_direction","NES Direction:",
                                  choices = c("All", "Positive", "Negative"), inline = TRUE),
                     downloadButton("download_table", "Download current table")
                     ),
                   mainPanel(
                     dataTableOutput("gsea_table")
                   )
                 )),
        tabPanel("NES Scatter Plot",
                 sidebarLayout(
                   sidebarPanel(
                     sliderInput("padj_cutoff_scatter", "Filter table by adjusted p-value (padj):",
                                 min = 0, max = 0.1, step = 0.01, value = 0.05)),
                   mainPanel(
                     plotOutput("nes_scatter")
                   )
                 ))
        )
      )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  # upload new fgsea dataset on submission
  load_fgsea_data <- reactive({
    req(input$submit_fgsea)
    req(input$fgsea_file)
    counts_df <-read.csv(input$fgsea_file$datapath)
    return(fgsea_res)
  })
  
  # render top pathways barplot
  output$pathways_barplot <- renderPlot({
    fgsea_res <- load_fgsea_data()
    top_pathways(fgsea_res, input$top_n)
  })
  
  # render filtered fgsea data table
  output$gsea_table <- renderDataTable({
    filtered_fgsea <- filter_fgsea(fgsea_res,input$padj_cutoff, input$nes_direction)
  })
  
  # render NES scatter plot
  output$nes_scatter <- renderPlot({
    plot_nes_scatter(fgsea_res,input$padj_cutoff_scatter)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
