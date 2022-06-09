net_to_ampl <-function(measurements,network)
  
{
  
  library(dplyr)
  library(stringr)
  library(tibble)
  load(measurements)
  load(network)
  
  measurements <- as.data.frame(TF_activities)
  network <- as.data.frame(network)
  
  network <- network %>% dplyr::mutate(Source = str_replace_all(as.character(Source), "\\.", "_"))
  network <- network %>% dplyr::mutate(Target = str_replace_all(as.character(Target), "\\.", "_"))
  source <- network[,1]
  sign <- network[,2]
  target <- network[,3]
  measurements <- add_column(measurements, Genes = rownames(measurements))
  
  measurements <- measurements %>% dplyr::mutate(Genes = str_replace_all(as.character(Genes), "\\.", "_"))
  Genes <- unique(unlist(measurements$Genes))
  
  set_j <- unique(unlist(c(source,target)))
  
  sink("carnival.dat")
  
  cat("param alpha := 5;")
  cat("\n")
  cat("param beta := 1;")
  cat("\n")
  cat("param p :=")
  cat("\n")
  cat("  TP53 1;")
  cat("\n")
  cat("set j := ")
  cat(set_j)
  cat(";")
  cat("\n")
  cat("set vv := ")
  cat(set_j)
  cat(";")
  cat("\n")
  cat("set v := ")
  cat(set_j)
  cat(";")
  cat("\n")
  
  cat("\n")
  
  cat("set IG:= ")
  cat("\n")
  
  
  for (i in 1: length(source)) {
    
    cat(paste0("(",source[i]),",",target[i],") \n")
    
  }
  cat(";")
  cat("\n")
  cat("param sigma := ")
  cat("\n")
  for (i in 1: length(source)) {
    
    cat(paste0(source[i])," ",target[i]," ",sign[i]," \n")
    
  }
  cat(";")
  
  
  
  cat("\n")
  cat("param m := ")
  cat("\n")
  for (i in 1: nrow(measurements)) {
    if (Genes[i] %in% c(source,target)) {
      cat(paste0(Genes[i])," ",measurements[i,1]," \n")
    }
  }
  cat(";")
  sink()
  
  
  
  return ;
}