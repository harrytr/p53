p53_mutations <- function() {
  
  
  library(ggplot2)
  library(dplyr)
  
  hotspots <- c("p.R175H","p.R248Q","p.R273H","p.R248W", "p.R273C", "p.R282W", "p.G245S")
  hotspots_tcga <- c("R175H","R248Q","R273H","R248W", "R273C", "R282W", "G245S")
  
  p53_CCLE <- as.data.frame(read.csv("p53_CCLE_all.csv"), header=TRUE)
  
  p53_TCGA <- as.data.frame(read.csv("p53_TCGA_all.csv"), header=TRUE)
  
  
  p53_TCGA  <-  p53_TCGA %>% dplyr::group_by(Variant_Classification, isDeleterious) 
  
  p53_TCGA <- p53_TCGA %>% dplyr::summarize(count = n())
  
  
  write.csv(p53_TCGA,"TCGA_mosaic.csv")
  
  
  p53_CCLE  <-  p53_CCLE %>% dplyr::group_by(Variant_Classification, isDeleterious) 
  
  p53_CCLE <- p53_CCLE %>% dplyr::summarize(count = n())
  
  
  write.csv(p53_CCLE,"CCLE_mosaic.csv")
  
  
  
  
  
  
}