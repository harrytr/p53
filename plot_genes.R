plot_genes <- function(df2) {

  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )
  df2 <- head(df2,100)
  
  
  gene <- "AGO2"
  print("Loading libraries required...")
  list.of.packages <- c("devtools","dplyr","ggplot2","ggrepel","ggpubr","viridis","tibble","stringr",
                        "corrplot","tidyverse","igraph","visNetwork", "data.table", "CARNIVAL",
                        "viper", "CellNOptR","edgeR", "OmnipathR", "stringi","openxlsx","samr",
                        "sna", "gplots","ggfortify","limma", "UpSetR")
  
  invisible(lapply(list.of.packages, library, character.only = TRUE))
  
  sp <- ggplot(data = df2,aes(x=factor(row.names(df2), levels=row.names(df2)), y = pvalues ))+
    #scale_fill_viridis_d( option = "D")+
    
    geom_point(size = 2,color = dplyr::case_when(df2$pvalues > 0.05 ~ "#FF0000",
                                                 df2$pvalues < 0.05 ~ "#00CC00",
                                                 TRUE ~ "#00CC00"), alpha = 0.8) +
    geom_hline(yintercept = 0.05, color = "#00CC00") +
    
   # geom_label_repel(aes(label=row.names(df2)),
    #                 box.padding   = 0.5,
     #                point.padding = 0.005, size = 2) +
    
    
    ylab(  c("P-value of Fisher exact test")  )  +
    xlab(  c("Hugo Symbol") ) +
    
    font("xylab",size=10)+
    font("xy",size=10)+
    font("xy.text", size = 10) +
    font("legend.text",size = 10) +
    theme(axis.text.x=element_text(size=5, angle=90,hjust=0.95,vjust=0.1)) +
    ggtitle(paste0("Genome-wide comparison of ", gene, " amplification versus MT/WT"))
  
  print(sp)
}