library(tidyverse)
library(ggplot2)
library('RColorBrewer')


#' function to load data
#' @param filename: file path to .csv file
#' @return dataframe: counts data with gene IDs column named "gene"
load_counts_data <- function(filename){
  counts_df <-read.csv(filename)
  colnames(counts_df)[1] <- "gene"
  return(counts_df)
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
    filter(gene_vars >= var_cutoff, gene_nonzero >= nonzero_filter) %>%
    select(gene)
  
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
  heatmap(heatmap_df, 
          col=brewer.pal(11, "RdBu"), 
          main = "Clustered Heatmap of Filtered Counts",
          cexRow=.5)
  # Plot a legend in bottom right part of heatmap
  legend(x = "bottomright", 
         legend = c("low", "medium", "high"),
         fill= colorRampPalette(brewer.pal(11, "RdBu"))(3)
         )
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
    select(c(pc1, pc2)) %>% # Extract cols for chosen principal components
    mutate(group = ifelse(grepl("^C", rownames(pca_data)), "control", "HD")) #create a group column based on the rownames
  
  # Create pca plot
  pca_plot <- ggplot(data = plotdata, aes(x = plotdata[,1], y = plotdata[,2], color = group)) +
    geom_point() +
    scale_color_discrete(name = "Sample Group") +
    labs(x=paste0("PC", pc1,": ", round(pca_var[pc1]*100,2),"% variance"),
         y=paste0("PC", pc2,": ", round(pca_var[pc2]*100,2),"% variance"))+
    theme_minimal()
  
  return(pca_plot)
}

  

