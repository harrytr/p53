radar_plot <- function() {
  
library(fmsb)

library(dplyr)

  
  
  
  df <- as.data.frame(read.csv("radar_plot_data_temp.csv"), header=TRUE)
  
  radarchart(df)
  
  
}