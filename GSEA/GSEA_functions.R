library('fgsea')
library(tidyverse)
library(ggplot2)

#' function to load deseq data
#' @param filename: file path to .csv file
#' @return dataframe: counts data with gene IDs column named "gene"
load_de_data <- function(filename){
  deseq_res <-read.csv(filename)
  colnames(deseq_res)[1] <- "ID"
  return(deseq_res)
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