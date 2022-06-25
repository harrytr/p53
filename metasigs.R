metasigs <- function(database) {
  
  library(UpSetR)
  keyword <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
  setwd(paste0("~","/",database))
  files <- list.files(getwd())
  print(files)
  all_sigs <- list()
  
  for (i in 1: length(files)) {
    print(files[i])
    df_primal <- as.data.frame(read.csv(files[i]), header = TRUE)
    
    sig <- as.list(unique(df_primal[2]))
    all_sigs[i] <- as.list(sig)
    names(all_sigs)[i] <- paste0(database,"_",files[i])

  }

  upset(fromList(all_sigs),order.by = "freq")
  
  
}