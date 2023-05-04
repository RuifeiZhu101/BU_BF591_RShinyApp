
#' function to load data
#' @param filename: file path to .csv file
#' @return dataframe: counts data with gene IDs column named "gene"
load_de_data <- function(filename){
  deseq_res <-read.csv(filename)
  colnames(deseq_res)[1] <- "ID"
  return(deseq_res)
}


#' Volcano plot
#'
#' @param dataf The loaded data frame.
#' @param x_name The column name to plot on the x-axis
#' @param y_name The column name to plot on the y-axis
#' @param slider A negative integer value representing the magnitude of
#' p-adjusted values to color. Most of our data will be between -1 and -300.
#' @param color1 One of the colors for the points.
#' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
#'
#' @return A ggplot object of a volcano plot

volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
  volcano<-ggplot(data = dataf,aes(x = !!sym(x_name), y=-log10(!!sym(y_name)))) +
      geom_point(aes(color = padj< 1*10^(slider)))+
      theme(legend.position = "bottom") +
      labs( color = str_glue('{y_name} 1 x 10^ {slider}'))
    return(volcano)
}

#' Draw and filter table
#'
#' @param dataf Data frame loaded by load_data()
#' @param slider Negative number, typically from the slider input.
#'
#' @return Data frame filtered to p-adjusted values that are less than 
#' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
#' displayed.
#' @details Not only does this function filter the data frame to 
#' rows that are above the slider magnitude, it should also change the format 
#' of the p-value columns to display more digits. This is so that it looks 
#' better when displayed on the web page. 

draw_table <- function(dataf, slider) {
  filtered_df <- dplyr::filter(dataf,padj<10^slider)
  #format the p-value columns to display more digit
  filtered_df$pvalue <- formatC(filtered_df$pvalue, format = "e", digits = 5)
  filtered_df$padj <- formatC(filtered_df$padj, format = "e", digits = 5)
  
  return(filtered_df)
}

