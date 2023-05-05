## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)
library(tidyverse)
source("counts_functions.R")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Counts Matrix Exploration"),
    p("Exploring and visualizing counts matrices can aid in 
    selecting gene count filtering strategies and understanding counts data structure.
      This component allows the user to choose different gene filtering thresholds
      and assess their effects using diagnostic plots of the counts matrix."),
 
    sidebarLayout(
        sidebarPanel(
          fileInput("counts_file",
                    "Load normalized counts matrix in CSV format:",
                    placeholder = "norm_counts_matrix.csv",
                    accept = ".csv"),
          #Slider to include genes with at least X percentile of variance
          sliderInput("var_cutoff",
                      "Slider to include genes with at least X percentile of variance:",
                      min = 0,
                      max = 100,
                      value = 50,
                      step = 1),
          #Slider to include genes with at least X non-zero samples
          sliderInput("nonzero_cutoff",
                      "Slider to include genes with at least X samples that are non-zero:",
                      min = 0,
                      max = 69, #total number of samples
                      value = 10,
                      step = 1),
          actionButton("submit_counts","Submit")
          ),

        mainPanel(
          #Tabs with summaries, plots, and heatmaps
          tabsetPanel(
            tabPanel("Summary",
                     tableOutput("summary_table")),
            tabPanel("Diagnostic Scatter Plots",
                     plotOutput("median_vs_variance"),
                     plotOutput("median_vs_nonzero")),
            tabPanel("Heatmap",
                     radioButtons("log_trans",
                                  "Enabling log-transforming counts for visualization:",
                                  choices = c("TRUE", "FALSE")),
                     plotOutput("heatmap")
                     ),
            tabPanel("PCA Plot",
                     numericInput("pc1",
                                  "PC:",
                                  value = 1,
                                  min = 1, 
                                  max = 69, 
                                  step = 1),
                     numericInput("pc2",
                                  "PC:", 
                                  value=2,
                                  min = 1,
                                  max = 69, 
                                  step = 1),
                     plotOutput("pca_plot"))
          )
        )
    )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) #increase the upload file limit to 30MB
  
  #' load counts data on submission
  #' rename the first column with ensembl gene IDs to "gene"
  load_counts_data <- reactive({
    req(input$submit_counts)
    req(input$counts_file)
    counts_df <-read.csv(input$counts_file$datapath)
    colnames(counts_df)[1] <- "gene"
    return(counts_df)
  })

  # render summary table
  output$summary_table <- renderTable({
    #load counts data
    counts_df <- load_counts_data()
    filtered_df <- filter_data(counts_df, input$var_cutoff, input$nonzero_cutoff)
    # call `create_summary_table()` function to create a summary table 
    create_summary_table(counts_df,filtered_df)
  })
  
 # render diagnostic scatter plot
  output$median_vs_variance <- renderPlot({plot_scatter(load_counts_data(),
                                                        input$var_cutoff,
                                                        input$nonzero_cutoff)[1]
    })
  output$median_vs_nonzero <- renderPlot({plot_scatter(load_counts_data(),
                                                        input$var_cutoff,
                                                        input$nonzero_cutoff)[2]
  })
  
  # render clustered heatmap
  output$heatmap <- renderPlot({plot_heatmap(filtered_df,
                                             input$log_trans)},
                               width = 600, height = 600)
  # render pca plot
  output$pca_plot <- renderPlot({plot_pca(filtered_df,
                                          input$pc1,
                                          input$pc2)
    })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
