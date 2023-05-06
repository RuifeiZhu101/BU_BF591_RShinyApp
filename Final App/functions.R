## Author: Ruifei Zhu
## BU BF591 final project
## Functions
library(tidyverse)
library(ggplot2)
library(dplyr)
library('fgsea')
library(gplots)
library(RColorBrewer)

#' Function to create a summary table for sample information
#'
#' @param data: sample information data
#' @param gmt (str): the path to the GMT file

#' @return table with sample information summary
sample_summary <- function(data){
  summary_data <- data.frame(
    Column_Name = names(data),
    Type = sapply(data, class),
    Mean_or_Distinct_Values <- sapply(data, function(x) {
      if (is.numeric(x) | is.integer(x)) {
        sprintf("%.2f (+/- %.2f)", mean(x,na.rm = TRUE), sd(x,na.rm = TRUE))
      } else if (is.factor(x) | is.character(x)) {
        paste(unique(x), collapse = ", ")
      } else {""}
    })
  )
  return(summary_data)
}


#' Function to create a histogram of continuous variables
#'
#' @param df: sample information data
#' @param x_var: variable to plot 
#' @param group_var: variable to group by

#' @return ggplot histogram of data counts
#' 
plot_sample_hist <- function(df,x_var,group_var){
  p <- ggplot(df, aes_string(x = x_var)) +
    geom_histogram(aes(fill = .data[[group_var]]), 
                   alpha = 0.5, 
                   position = "identity", 
                   bins = 30,
                   na.rm = TRUE) +
    theme_bw() +
    labs(x = x_var, y = "Count", fill = group_var) 
  p
}


#' function to filter data based on percentile of variance and number of zeros
#' 
#' @param counts_df: normalized counts data
#' @param var_filter: minimum percentile of variance 
#' @param nonzero_filter: minimum number of zero samples for a gene
#' 
#' @return data.frame: a filtered dataframe with a 'gene' column followed by
#' sample names as column names.
#' 
filter_data <- function(counts_df, var_filter, nonzero_filter){
  df <- counts_df %>%
    pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
    group_by(gene) %>% 
    summarise(gene_vars = var(counts), gene_nonzero = sum(counts != 0)) %>% 
    ungroup() 
  
  # Calculate variance cutoff based on var_filter percentile
  var_cutoff <- quantile(df$gene_vars, var_filter/100, na.rm = TRUE)
  
  # Filter genes based on variance and nonzero sample cutoffs
  filtered_gene <- df %>%
    dplyr::filter(gene_vars >= var_cutoff, gene_nonzero >= nonzero_filter) %>%
    dplyr::select(gene)
  
  # Filter original data based on selected genes
  filtered_df<- counts_df %>% filter(gene %in% filtered_gene$gene)
  
  return(filtered_df)
}

#' function to filter data based on percentile of variance and number of zeros
#' 
#' @param counts_df: normalized counts data
#' @param filtered_df: a filtered dataframe according to a percentile of variance and number of zeros
#' 
#' @return data.frame: a dataframe with the summary of sample and gene counts statistics
#'
create_summary_table <- function(counts_df, filtered_df){
  # calculate summary statistics
  total_genes <- nrow(counts_df)
  total_samples <- ncol(counts_df) - 1
  filtered_genes <- nrow(filtered_df) 
  nopass_genes <- total_genes - filtered_genes
  
  summary_df <- data.frame(
    Total_Samples = total_samples,
    Total_Genes = total_genes,
    Filtered_Genes = filtered_genes,
    Nopass_Genes = nopass_genes,
    Percent_Filtered = paste0(round(filtered_genes / total_genes * 100,2),"%"),
    Percent_Nopass = paste0(round(nopass_genes / total_genes * 100,2),"%")
  )
  summary_df <-data.frame(lapply(summary_df, as.character), stringsAsFactors=FALSE)
  colnames(summary_df) <- c("Total Samples", "Total Genes","Filtered Genes", "Genes not passing filter", "% Filtered genes","% Genes not passing filter")
  return(summary_df)
}

#' function to create diagnostic scatter plot  based on percentile of variance and number of zeros
#' 
#' @param counts_df: normalized counts data(data.frame,with gene ID column names "gene")
#' @param var_filter: minimum percentile of variance 
#' @param nonzero_filter: minimum number of zero samples for a gene
#' 
#' @return ggplot: two scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter: 
#' p_var_median: Variance vs. median count (log10), 
#' p_median_zeros: median count vs number of zeros
#' 
plot_scatter <- function(counts_df, var_filter, nonzero_filter) {
  filtered_df <- filter_data(counts_df, var_filter, nonzero_filter)
  p_var_median<- counts_df %>% 
    pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
    group_by(gene) %>% 
    summarise(gene_median = median(counts), gene_var = var(counts))%>% 
    ungroup() %>% 
    mutate(pass_filter = ifelse(gene %in% filtered_df$gene, "Yes", "No")) %>% 
    ggplot(., aes(x = log10(gene_median + 1), y = log10(gene_var), colour = pass_filter)) +
    # add a pseudocount to avoid log transformation issues with 0 counts
    geom_point() +  
    scale_colour_manual(values = c("No" = "grey", "Yes" = "black")) + 
    labs(title = "Variance vs. Median Counts",
         x = "Gene count median(log10 scale)",
         y = "Gene count variance(log10 scale)") +
    theme_bw()
  
  
  p_median_zeros<- counts_df %>% 
    pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
    group_by(gene) %>% 
    summarise(gene_median = median(counts), gene_zeros = sum(counts== 0))%>% 
    ungroup() %>% 
    mutate(pass_filter = ifelse(gene %in% filtered_df$gene, "Yes", "No")) %>% 
    ggplot(., aes(x = log10(gene_median + 1), y = gene_zeros, colour = pass_filter)) +
    geom_point() +  
    scale_colour_manual(values = c("No" = "grey", "Yes" = "black")) + 
    labs(title = "number of zeros vs. Median Counts  ",
         x = "Gene count median(log10 scale)",
         y = "Number of samples with zero count") +
    theme_bw()
  return(list(p_var_median,p_median_zeros))
}

#' A function to create a clustered heatmap of counts remaining after filtering
#'
#' @param filtered_df: a filtered dataframe according to a percentile of variance and number of zeros
#' @param log_trans: boolean: whether to transform counts to log10 values. 
#'
#' @return A clustered heatmap displaying the counts for the differentially expressed genes

plot_heatmap <- function(filtered_df,log_trans) {
  # remove the "gene" column and transpose the dataframe to 
  # display genes on the y-axis and the samples on the x-axis of the heatmap 
  heatmap_df <- as.matrix(filtered_df[, -1])
  
  # transform counts to log10 values if log_trans is TRUE
  if(log_trans) {
    heatmap_df <- log10(heatmap_df + 1)
  } else {
    heatmap_df <- heatmap_df
  }
  # create heatmap with color bar
  heatmap.2(heatmap_df, 
          col=brewer.pal(11, "YlOrRd"), 
          main = "Clustered Heatmap of Filtered Counts")
}

#' A function to create a scatter plot of principal component analysis projections 
#' allow the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2)
#'
#' @param filtered_df: a filtered dataframe according to a percentile of variance and number of zeros
#' @param pc1: integer, first chosen principal component to plot
#' @param pc2: integer, second chosen principal component to plot
#' @return A scatter plot of principal component analysis

plot_pca <- function(filtered_df,pc1,pc2) {
  # Perform PCA analysis on filtered data frame
  pca_results <- prcomp(scale(t(filtered_df[,-1])), center=FALSE, scale=FALSE)
  pca_data <- as.data.frame(pca_results$x)
  pca_var <- summary(pca_results)$importance[2,] #variance of all PCs
  # create plot data 
  plotdata <- pca_data %>% 
    dplyr::select(c(pc1, pc2)) %>% # Extract cols for chosen principal components
    dplyr::mutate(group = ifelse(grepl("^C", rownames(pca_data)), "control", "HD")) #create a group column based on the rownames
  
  # Create pca plot
  pca_plot <- ggplot(data = plotdata, aes(x = plotdata[,1], y = plotdata[,2], color = group)) +
    geom_point() +
    scale_color_discrete(name = "Sample Group") +
    labs(x=paste0("PC", pc1,": ", round(pca_var[pc1]*100,2),"% variance"),
         y=paste0("PC", pc2,": ", round(pca_var[pc2]*100,2),"% variance"))+
    theme_minimal()
  
  return(pca_plot)
}


#' Volcano plot
#'
#' @param dataf The loaded data frame.
#' @param x_name The column name to plot on the x-axis
#' @param y_name The column name to plot on the y-axis
#' @param slider A negative integer value representing the magnitude of
#' p-adjusted values to color. Most of our data will be between -1 and -300.
#' @param color1 One of the colors for the points below the
#' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
#'
#' @return A ggplot object of a volcano plot

volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
  #clean data
  dataf <- dataf %>%
    drop_na(!!sym(x_name)) %>%
    drop_na(!!sym(y_name))
  #create plot
  p <- dataf %>%
    ggplot(aes(!!sym(x_name), -log10(!!sym(y_name)), 
               colour = log10(!!sym(y_name)) < slider)) +
    geom_point(stat="identity")+
    scale_colour_manual(paste0("P-adj < 1e",slider),values = c("TRUE" = color2, "FALSE" = color1)) +
    xlab(x_name)+
    ylab(y_name)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          legend.position="bottom")
  p
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


#' Function to run fgsea on DESeq2 results
#'
#' @param deseq_res: differential expression analysis results from DESeq2
#' @param gmt (str): the path to the GMT file

#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(deseq_res, 'c2.cp.v7.5.1.symbols.gmt')
run_gsea <- function(deseq_res, gmt) {
  # creating  a named vector [ranked genes by stat]
  res <- deseq_res[order(-deseq_res$log2FoldChange),] 
  gene_list <- res$log2FoldChange
  names(gene_list) <- res$symbol
  
  # get C2 gene sets
  C2_gene_sets <- fgsea::gmtPathways(gmt)
  
  # run fgsea
  fgsea_res <- fgsea(pathways = C2_gene_sets,
                     stats = gene_list,
                     minSize = 15,
                     maxSize = 500)
  
  return(fgsea_res)
}

#' Function to plot top n positive NES and top n negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  # select the top positive and negative pathways
  top_positive <- fgsea_results %>%
    arrange(desc(NES)) %>%
    slice_head(n = num_paths) %>%
    filter(NES > 0)
  
  top_negative <- fgsea_results %>%
    arrange(NES) %>%
    slice_head(n = num_paths) %>%
    filter(NES < 0)
  
  # combine the data into one dataframe
  top_paths <- rbind(top_positive, top_negative)
  # create the ggplot
  bar_plt<-ggplot(top_paths, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "blue"))+ 
    coord_flip() +
    labs(x = "", y = "Normalized Enrichment Score (NES)") +
    theme(axis.text = element_text(size = 1))+
    theme_classic()
  return(bar_plt)
}

#' Function to filter table by adjusted p-value
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param padj_cutoff(num): the threshold of adjusted p-value to filter out data with padj >= padj_cutoff
#' @param nes_direc(str): 'All','Positive','Negative',to select all, positive or negative NES pathways
#' @return filtered fgsea results data table 
#' @export
#'
#' @examples filtered_fgsea <- filter_fgsea(fgsea_results, 0.05, "All" )
filter_fgsea <- function(fgsea_results, padj_cutoff, nes_dirct){
  # filter by adjusted p-value cutoff
  fgsea_res_filtered <- fgsea_results %>% 
    filter(padj <= padj_cutoff)
  # filter by NES direction
  if (nes_dirct == "pos") {
    fgsea_res_filtered <- fgsea_res_filtered %>% 
      filter(NES > 0)
  } else if (nes_dirct == "neg") {
    fgsea_res_filtered <- fgsea_res_filtered %>% 
      filter(NES < 0)
  }
  return(fgsea_res_filtered )
}

#' Function to creat scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, 
#' with gene sets below threshold in grey color
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param padj_cutoff(num): the threshold of adjusted p-value to filter out data with padj >= padj_cutoff
#' 
#' @return ggplot: scatter plot  of NES vs. -log10 adjusted p-value

plot_nes_scatter  <- function(fgsea_results, padj_cutoff){
  # filter by adjusted p-value cutoff
  fgsea_res_filtered <- fgsea_results %>% 
    filter(padj <= padj_cutoff)
  
  # add a color column for plot data
  pltdata <- fgsea_results%>%mutate(color = ifelse(padj <= padj_cutoff, "Yes", "No"))
  # create scatter plot
  p <- ggplot(pltdata, aes(x = NES, y = -log10(padj), color = color)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(name = paste0("padj <", padj_cutoff),
                       values = c("grey", "red"),
                       labels = c("No", "Yes")) +
    theme_minimal() +
    labs(x = "Normalized Enrichment Score (NES)", y = "-log10(Padj)")
  
  return(p)
}




