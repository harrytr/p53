preprocess_hypoxia_1 <- function() {
  
  
  
  library(dplyr)
  library(data.table)
  library(tibble)
  
  anno <- as.data.frame(read.csv("sample_info1.csv", header=TRUE))
  
  RNA <- as.data.frame(read.csv("genecounts1.csv", header=TRUE))
  
  
  anno <- add_column(anno,new_id = "")
  #anno <- anno %>% unite("new_id","Filename","Sample",	"Cell_Line"	,"Treatment",	"Media"	,"Glucose"	,"Glutamaine"	,"Time", sep = "_", remove = FALSE)
  
  anno <- anno %>% unite("new_id","Filename",	"Sample",	"Count"	,"Uniquely_Mapped_Reads","Cell_Line",	"Treatment",	"Media",	"Batch")
  
  anno <- anno %>% select("new_id")
  
  
  #colnames(anno)[2] <- "Tumor_Sample_ID"
  colnames(RNA)[2:length(colnames(RNA))] <- anno$new_id 
  
  write.csv(RNA,"RNA_hypoxia.csv")
  
}