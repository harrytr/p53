bc_signatures_single_all <- function(database) {
  library(UpSetR)
  all_sigs <- list()
  setwd("C:/Users/Harry/OneDrive - University of Greenwich/Desktop/RESULTS/CENTRALITY_PREPROCESS")
  setwd(paste0(getwd(),"/",database))
  
  
  keyword <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense","Nonsense","Splice_Site","_True","_False")
  for (i in 1:length(keyword)) {
    print(i)
    sig<- bc_signatures_single(database,i)
    all_sigs[i] <- (sig)
    names(all_sigs)[i] <- paste0(database,"_",keyword[i])
    
    
  }
  upset(fromList(all_sigs),nsets = length(keyword),order.by = "degree")
}