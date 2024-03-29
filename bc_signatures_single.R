bc_signatures_single <- function(database,single) {
  
  
  keywords <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
  keyword <- c(keywords[single])


  #single <- 1
  library(dplyr)
  library(tidyverse)
  library(UpSetR)
  
  setwd("C:/Users/Harry/OneDrive - University of Greenwich/Desktop/RESULTS/CENTRALITY_PREPROCESS")
  setwd(paste0(getwd(),"/",database))
  files <- list.files(getwd())
  df2  <- matrix(data = , nrow = 1, ncol = 3) #
  df_final  <- matrix(data = , nrow = 0, ncol = 3) #
  colnames(df2) <- c("Cancer_type", "Mutation_Type/Deleterious", "Signature")
  colnames(df_final) <- c("Cancer_type", "Mutation_Type/Deleterious", "Signature")
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
  #row_no <- 0
  size <- length(files)
  
  for (i in 1:(size-1)) {
    df2  <- matrix(data = , nrow = length(keyword), ncol = 3) #
    df <- matrix(data = , nrow = , ncol = 1) #
    df_primal <- as.data.frame(read.csv(files_bkp[i]), header = TRUE)
    for (j in 1:length(keyword)) {
      df2[j,1] <- files[i]
      df2[j,2] <- keyword[j]
      print(keyword[j])
      df <- df_primal  %>%  dplyr:: filter(str_detect(df_primal$top_bc_network ,keyword[j])==TRUE)
      print(df)
      if (nrow(df)>0) {
        temp1 <- sapply(str_split(df$top_bc_vertex,"_"),function(x) ((x)))
        temp <- sapply(unlist(temp1),function(x) ((x)))
        signature <- unique(names(temp))
        signature <- paste(signature, collapse = ",")
        df2[j,3] <- signature
      }
      else{
        print(files[i])
        signature <- ""
      }
      #g <- upset(df2, order.by = "freq")
      #print(g)
    }
    print(files[i])
    print("-----------------------")
    df_final <- rbind(df_final,df2)
    
  }
  setwd(paste0("C:/Users/Harry/OneDrive - University of Greenwich/Desktop/RESULTS/CENTRALITY_PREPROCESS/meta-sigs/",database))
  df_final <- as.data.frame(df_final)
  df_final <- df_final %>% dplyr::filter(Signature != "NA")
  write.csv(df_final,paste0(database,"_",keywords[single],"_signatures.csv"))
  
    
  signature <- (paste(df_final[,3], collapse = ","))
  temp1 <- sapply(str_split(signature,","),function(x) ((x)))
  signature <- as.list(temp1)
  print(length(signature))
  signature <- unique(signature)
  print(length(signature))
  signature <- list(signature)
 
  write.csv(signature,paste0(keywords[single],".csv"))
 
  return(signature)
}