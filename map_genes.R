map_genes <- function() {
  
  library(dplyr)
  
  #df1 <- as.data.frame(read.csv("genome_results_breast.csv", header=TRUE)) #
  df1 <- as.data.frame(read.csv("CCLE_MT_VS_WT_adjpvalues_breast.csv", header=TRUE))
  colnames(df1)[1] <- "gene"
  
  #df2 <- as.data.frame(read.csv("genome_results_metabric.csv", header=TRUE)) #
  df2 <- as.data.frame(read.csv("METABRIC_MT_VS_WT_adjpvalues.csv", header=TRUE))
  
  colnames(df2)[1] <- "gene"
  #df3 <- as.data.frame(read.csv("genome_results_brca.csv", header=TRUE)) # 
  df3 <- as.data.frame(read.csv("TCGA_BRCA_MT_VS_WT_adjpvalues.csv", header=TRUE))
  
  colnames(df3)[1] <- "gene"
  
  #signature <- df1 %>% dplyr::filter(df1$pvalues < 0.05)
  signature <- df1$gene
  #signature <- signature[,1]
  print(signature)
  
  
  sub_metabric <- df2 %>% dplyr::filter(gene %in% signature)

  sub_tcga <- df3 %>% dplyr::filter(gene %in% signature)
  
  write.csv(sub_metabric,"Sub_metabric.csv")
  write.csv(sub_tcga,"Sub_tcga.csv")
  
}