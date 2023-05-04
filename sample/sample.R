## Author: Ruifei Zhu
## BU BF591
## Final Project - RShiny application that features multiple bioinformatics processes implemented in R

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(DT)



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Sample Information Exploration"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          fileInput("sampleinfo_file",
                    "Load sample information matrix in CSV format:", 
                    placeholder = "sample_metadata.csv", 
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
                                      choices = ""),
                         actionButton(inputId = "plot_choice",
                                      label = "Plot"),
                         icon = icon("redo")
                         ),
                       mainPanel(plotOutput("histogram"))
                       )
                     )
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #' load_Data on submission
  load_data <- reactive({
    req(input$submit_sampleinfo)
    req(input$sampleinfo_file)
    sample_info <-read.csv(input$sampleinfo_file$datapath)
    return(sample_info)
  })
  
  # render summary table
  output$summary_table <- renderDataTable({
    data <- load_data()
    summary_data <- data.frame(
      Column_Name = names(data),
      Type = sapply(data, class),
      Mean_or_Distinct_Values <- sapply(data, function(x) {
        if (is.numeric(x) | is.integer(x)) {
          sprintf("%.2f (+/- %.2f)", mean(x,na.rm = TRUE), sd(x,na.rm = TRUE))
        } else if (is.factor(x) | is.character(x)) {
          paste(unique(x), collapse = ", ")
        } else {
          ""
        }
      })
    )
    DT::datatable(summary_data,
              options = list(pageLength = 10, lengthMenu = c(5, 10, 15)),
              colnames = c('Column Name', 'Type', 'Mean (sd) or Distinct Values'),
              rownames = FALSE)
  })
  # render data table
  output$data_table <- renderDataTable({
    df <- load_data()
    DT::datatable(df,
                  options = list(pageLength = 10, 
                                 lengthMenu = c(5, 10, 15, 20),
                                 scrollX = TRUE) # limit the horizontal growth of the page
                                                 #scroll to see all the columns
                  )
  })
  
  # update plot radio buttons on data input
  observe({
    df <- load_data()
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
    df <- load_data()
    x_var <- input$plot_col
    group_var <- input$group_col
    df[[group_var]] <- factor(df[[group_var]])
    p <- ggplot(df, aes_string(x = x_var)) +
      geom_histogram(aes(fill = .data[[group_var]]), 
                     alpha = 0.5, 
                     position = "identity", 
                     bins = 30,
                     na.rm = TRUE) +
      theme_bw() +
      labs(x = x_var, y = "Count", fill = group_var) 
    print(p)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
