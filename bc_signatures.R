bc_signatures <- function(database) {
  
  
    keyword <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
    library(dplyr)
    library(tidyverse)

    setwd(paste0("~","/",database))
    files <- list.files(getwd())
    df2  <- matrix(data = , nrow = length(keyword), ncol = 3) #
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
      }
      print(files[i])
      print("-----------------------")
      df_final <- rbind(df_final,df2)
     
    }

    write.csv(df_final,paste0(database,"_signatures.csv"))

    return()
  }