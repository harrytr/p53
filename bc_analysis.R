bc_analysis <- function(database) {
  
  keyword <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
  library(dplyr)
  library(tidyverse)
  

  
  setwd(paste0("/Users/ctrianta/Documents/shinyapp_2022/centrality/",database))
  files <- list.files(getwd())
  df2 <- matrix(data = , nrow = length(files), ncol = length(keyword)) #
  colnames(df2) <- keyword
  files_bkp <- files
  
  
  
  files <- sapply(strsplit(as.character(files), split='.csv.csv', fixed=TRUE),function(x) (x[1]))
  if (database == "CCLE") {
    
    files <- sapply(strsplit(as.character(files), split='_', fixed=TRUE),function(x) (x[3]))
  }
    
    else
    {
  files <- sapply(strsplit(as.character(files), split='_', fixed=TRUE),function(x) (paste0(x[3],"_",x[4])))
  
    }
  files <- sapply(strsplit(as.character(files), split='.', fixed=TRUE),function(x) (x[1]))
  disease <- c()

  for (i in 1:length(files)) {
    print(files[i])
    for (j in 1:length(keyword)) {
      temp_name <- files[i]
      
      df <- as.data.frame(read.csv(files_bkp[i]), header = TRUE)
      colnames(df)[1] <- "central_gene"
      
      df <- df  %>%  dplyr:: filter(str_detect(df$top_bc_network ,keyword[j])==TRUE)

      central_gene <- which.max(table(df$central_gene))
      central_gene <- names(central_gene)
      print(central_gene)
      if (is.null(central_gene)) {central_gene <- "-"}
      df2[i,j] <- central_gene
    }
    disease <- c(disease,temp_name)
  }
  df_final <- cbind(disease,df2)
  #print(df_final)
  colnames(df_final)[1] <- paste0(database)
 # write.csv(df_final,paste0("/Users/ctrianta/Desktop/centrality/centrality_genes_",database,".csv"))
  unique_genes <- c()
  
  for (i in 2: ncol(df_final)) {
    unique_genes <- c(unique_genes, unique(df_final[,i]))
    
  }
  print("List of genes:")
  print(unique(unique_genes))
  return()
}