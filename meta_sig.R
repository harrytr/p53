meta_sig <- function()
{
  setwd("/rds/general/project/piel-chadeau_sickle_cell_disease/live/Copy of Original data/Phase1/sd2anon/csv/AsthmaUBIOPRED/CCLE_TCGA_SIGS/meta-sigs")
  df <- as.data.frame(read.csv("CCLE_Frame_Shift_Del_signatures.csv"), header = TRUE)
  print(df$Signature)
  df <- df %>% dplyr::filter(!is.na(df$Signature))

  #temp1 <- sapply(str_split(df$signature,","),function(x) ((x)))
  #temp <- sapply(unlist(temp1),function(x) ((x)))
  #signature <- unique(names(temp))

  signature <- unique(paste(df$Signature, collapse = ","))

  return(signature)
  
  
  
}