## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)
library(tidyverse)
library(dplyr)
source("functions.R")



# User Interface
ui <- fluidPage(
    # Application title
    titlePanel("BF591 Final Project - R Shiny Application for RNA-seq data exploration"),
    p("A demo dataset to help you understand the use of this application can be downloaded from the data directory of this app's repository."),
    # Create a tabset that includes four tabs: Sample, Counts, DE, GSEA
    tabsetPanel(
     ##-------Sample UI----
      tabPanel("Samples",
               h3("Sample Information Exploration"),
               p("The distinct values and distributions of sample information are important to understand before conducting analysis of corresponding sample data. 
                 This component allows the user to load and examine a sample information matrix."),
               p("Inputs: Sample information matrix in CSV format"),
               sidebarLayout(
                 sidebarPanel(
                   # Sample information Input: .csv file
                   fileInput("sampleinfo_file",
                             "Load sample information matrix in CSV format:", 
                             placeholder = "sample_info.csv", 
                             accept = ".csv"),
                   # Submit button
                   actionButton(inputId = "submit_sampleinfo",label = "Submit")
                   ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary",dataTableOutput("info_summary_table")),
                     tabPanel("Table", dataTableOutput("data_table")),
                     tabPanel("Plot",
                              # add a side bar to let user choose the column to plot
                              # and a column to group by
                              sidebarLayout(
                                sidebarPanel(
                                  radioButtons(inputId = "plot_col",
                                               label = "Select a varible to plot:",
                                               choices = ""),
                                                # to be gotten from the uploaded data file
                                  radioButtons(inputId = "group_col",
                                               label = "Select a group by condition:",
                                               choices = "")
                                  ),
                                mainPanel(plotOutput("histogram"))
                                 )
                              )) 
                   ),
                 ),
               ),
      
      ##------Counts UI---
        tabPanel("Counts",
                 h3("Counts Matrix Exploration"),
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
                                tableOutput("counts_summary_table")),
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
              ),
      ##----------DE UI---
        tabPanel("DE",
                 h3("Differential Expression"),
                 p("This component allows the user to load and explore a differential expression dataset."),
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
                     
                     colourInput("base","Base point color","#22577A"),
                     
                     colourInput("highlight","Highlight point color","#FFCF56"),
                     
                     sliderInput("slider",
                                 p("Select the magnitude of the p adjusted coloring:"),
                                 min = -30, max = 0, value = -10),
                     actionButton("submit_de", "Plot",icon = icon("redo"))
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
        ),
      ##-------- GSEA UI---
        tabPanel("GSEA",
                 h3("Gene Set Enrichment Analysis with fgsea"),
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
    )
  )

server <- function(input, output, session) {
  
  ##-------Sample Server---
  #' load metadata on submission
  load_metadata <- reactive({
    req(input$submit_sampleinfo)
    req(input$sampleinfo_file)
    sample_info <-read.csv(input$sampleinfo_file$datapath)
    return(sample_info)
  })
    
    # render sample information summary table
    output$info_summary_table <- renderDataTable({
      metadata <- load_metadata()
      summary_table <- sample_summary(metadata)
      DT::datatable(summary_table,
                options = list(pageLength = 10, lengthMenu = c(5, 10, 15)),
                colnames = c('Column Name', 'Type', 'Mean (sd) or Distinct Values'),
                rownames = FALSE)
    })
    
    # render data table
    output$data_table <- renderDataTable({
      metadata <- load_metadata()
      DT::datatable(metadata,
                    options = list(pageLength = 10, 
                                   lengthMenu = c(5, 10, 15, 20),
                                   scrollX = TRUE)# limit the horizontal growth of the page
                    #scroll to see all the columns
      )
    })
    
    # update plot radio buttons on data input
    observe({
      df <- load_metadata()
      continuous_cols <- sapply(df, function(x) is.numeric(x) | is.integer(x))
      updateRadioButtons(session, "plot_col", 
                         choices = names(df)[continuous_cols],
                         selected = names(df)[which(continuous_cols)[1]])
      updateRadioButtons(session, "group_col", 
                         choices = names(df),
                         selected = names(df)[1])
    })
    
    # render histogram plot
    output$histogram <- renderPlot({
      df <- load_metadata()
      x_var <- input$plot_col
      group_var <- input$group_col
      df[[group_var]] <- factor(df[[group_var]])
      plot_sample_hist(df,x_var ,group_var )
    })
    ##------------Counts server---
    
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
    output$counts_summary_table <- renderTable({
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
    output$heatmap <- renderPlot({
      counts_df <- load_counts_data()
      filtered_df <- filter_data(counts_df, input$var_cutoff, input$nonzero_cutoff)
      plot_heatmap(filtered_df,input$log_trans)},
      width = 600, height = 600)
    # render pca plot
    output$pca_plot <- renderPlot({
      counts_df <- load_counts_data()
      filtered_df <- filter_data(counts_df, input$var_cutoff, input$nonzero_cutoff)
      plot_pca(filtered_df,input$pc1,input$pc2)
    })
    
    #------------DE server---
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
      deseq_res <- load_de_data()
      draw_table(deseq_res,input$slider)
    })
    ##-----------GSEA server---
    # Use the reactive value to load the DE data and run fgsea
    fgsea_file <- reactive({
      deseq_res <- load_de_data()
      if (input$submit_fgsea && !is.null(input$fgsea_file)) {
        # Use uploaded FGSEA file
        fgsea_res <- read.csv(input$fgsea_file$datapath)
      } else {
        # Run FGSEA using DE seq results file
        fgsea_res <- run_gsea(deseq_res, 'c2.cp.v7.5.1.symbols.gmt')
      }
      return(fgsea_res)
    })
    

    # render top pathways barplot
    output$pathways_barplot <- renderPlot({
      fgsea_res <- fgsea_file()
      top_pathways(fgsea_res, input$top_n)
    })
    
    # render filtered fgsea data table
    output$gsea_table <- renderDataTable({
      fgsea_res <-  fgsea_file()
      filtered_fgsea <- filter_fgsea(fgsea_res,input$padj_cutoff, input$nes_direction)
    })
    
    # download filtered fgsea data table
    output$download_table <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        fgsea_res <- fgsea_file()
        filtered_fgsea <- filter_fgsea(fgsea_res, input$padj_cutoff, input$nes_direction)
        write.csv(filtered_fgsea, file)
      }
    )
    
    # render NES scatter plot
    output$nes_scatter <- renderPlot({
      fgsea_res <- fgsea_file()
      plot_nes_scatter(fgsea_res,input$padj_cutoff_scatter)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
