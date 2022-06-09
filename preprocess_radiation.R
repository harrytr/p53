preprocess_radiation <- function() {
  
  
  
  library(dplyr)
  library(data.table)
  
  anno <- as.data.frame(read.csv("H460_annotation.csv", header=TRUE))
  
  RNA <- as.data.frame(read.csv("H460_logcounts.csv", header=TRUE))
  
  
  anno <- add_column(anno,new_id = "")
  anno <- anno %>% unite("new_id", "ID","CELL.LINE",	"TIME",	"CLASS",	"nBATCH",	"SUB.LINE", sep = "_", remove = FALSE)
  anno <- anno %>% select("new_id", "ID")
  
  
  colnames(anno)[2] <- "Tumor_Sample_ID"
  colnames(RNA)[2:length(colnames(RNA))] <- anno$new_id 
  
  write.csv(RNA,"RNA.csv")
  
}