## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)



# User Interface
ui <- fluidPage(
    # Application title
    titlePanel("BF591 Final Project - R Shiny Application for RNA-seq data exploration"),
    p("A demo dataset to help you understand the use of this application can be downloaded from the data directory of this app's repository."),
    # Create a tabset that includes four tabs: Sample, Counts, DE, GSEA
    tabsetPanel(
     #<----------------------------Sample tab-------------------------------->
      tabPanel("Samples",
               h3("Sample Information Exploration"),
               p("The distinct values and distributions of sample information are important to understand before conducting analysis of corresponding sample data. 
                 This component allows the user to load and examine a sample information matrix."),
               p("Inputs: Sample information matrix in CSV format"),
               sidebarLayout(
                 sidebarPanel(
                   # Sample information Input: .csv file
                   fileInput("sampleinfo_file",
                             "Load sample information matrix in CSV format", 
                             placeholder = "sample_info.csv", 
                             accept = ".csv"),
                   # Submit button
                   actionButton(inputId = "submit_sampleinfo",label = "Submit")
                   ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary",dataTableOutput("summary_table")),
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
      #<----------------------------Counts tab------------------------------>
        tabPanel("Counts"),
      #<----------------------------DE tab------------------------------>
        tabPanel("DE"),
      #<---------------------------- One more analysis tab------------------------------>
        tabPanel("")
      ),
      
    )


# Define server logic required to draw a histogram
server <- function(input, output) {
  load_metadata <- reactive({
    input$submit_sampleinfo
    # Allow changes on submission
    isolate({
      sample_info <- read.csv(input$sampleinfo_file$datapath)
      #sample_info <- read.csv("data/sample_information.csv")
      return(sample_info)
    })
    
  })

    
}

# Run the application 
shinyApp(ui = ui, server = server)
