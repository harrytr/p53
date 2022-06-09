 p_adj <- function(csv) {
   
 library(tibble)
 df <- as.data.frame(read.csv(csv, header=TRUE))
   
   
 df <- add_column(df, adjusted_pvalue = 0)
   
 n = length(df$pvalues <= 0.05)
 df$adjusted_pvalue <- p.adjust(df$pvalues, method = "BH", n = n)
  
 
 write.csv(df,csv) 
 }

