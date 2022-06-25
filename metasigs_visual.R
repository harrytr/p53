metasigs_visual <- function(database) {
  
  keyword <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
  setwd(paste0("~","/",database))
  files <- list.files(getwd())
  
  for (i in 1: length(files)) {
    print(files[i])
    df_primal <- as.data.frame(read.csv(files[i]), header = TRUE)
    
    
    
    temp1 <- sapply(str_split(df_primal[2],","),function(x) ((x)))
    print(nrow(temp1))
    signature <- unique(temp1)
    
    print(length(signature))
    write.csv(temp1,paste0("set_",files[i]))
  }

  
  
  
  
  
}