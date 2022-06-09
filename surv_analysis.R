

surv_analysis <-function() {
  
  # using a cox proportional hazards model to analyse TCGA surv data
  library("survival")
  library("survminer")
  library("ggfortify")
  library("dplyr")
  library("tidyverse")
  splots <- list()
  
  keywords <- "R175H|R248Q|R273H|R248W|R273C|R282W|G245S"
  # keywords <- "_True"  #del
 # keywords <- "Missense"
  #keywords <- "WT"
  
  #feature <- "Missense_Del"
  feature <- "Hotspot-Deleterious"
  #inv_feature <- paste0("Non-", feature)
  inv_feature <- "Non-Hotspot & Deleterious"
  
  df2 <- as.data.frame(read.csv("opt_surv.csv", header=TRUE))
  
  # df2 <- df2  %>%  dplyr:: filter(OS.time <= 3000)
   df2 <- df2[-(grep("WT", df2$DEL_function)),]
  
  #  df2 <- df2  %>%  dplyr:: filter(type %in% "BRCA")
 df2 <- df2  %>%  dplyr:: filter(str_detect(df2$DEL_function,"_True")==TRUE)
  
  
  df2 <- mutate(df2 , DEL_function = ifelse(str_detect(df2$DEL_function,keywords)==TRUE,1,0))
  
  df2<-as.data.frame(df2[!(df2$OS.time=="#N/A"),])
  df2$OS.time <- as.numeric(df2$OS.time)
  df2$OS <- as.numeric(df2$OS)
  df2$DEL_function <- as.numeric(df2$DEL_function)
  df2$age_at_initial_pathologic_diagnosis <- as.numeric(df2$age_at_initial_pathologic_diagnosis)
  
  
  print(summary(coxph(Surv(OS.time, OS ) ~ DEL_function + strata(type), data=df2)))
  #stop()
  fit <- survfit(Surv(OS.time, OS ) ~ DEL_function, data = df2)
  
  gg <- ggsurvplot(fit,
                   font.main = 25,
                   font.x = 20,
                   font.y = 22,
                   font.tickslab = 15,
                   pval = TRUE,
                   break.time.by = 365,
                   censor.shape = "*", 
                   title = "Survival Curve", 
                   data = df2,
                   legend.labs= c(inv_feature,feature),
                   xlab = "Time (days)",
                   conf.int = TRUE,
                   risk.table = TRUE, 
                   risk.table.col = "strata",
                   
                   palette = c("#E7B800", "#2E9FDF"),
                   
                   risk.table.height = .3,
                   xlim = c(0,3650),
                   font.legend = c(25,"bold"),
                   #ncensor.plot.height = 0.25,
                   #conf.int.style = "step",  # customize style of confidence intervals
                   surv.median.line = "hv",  # add the median survival pointer.
  )
  
  print(gg)
  # ggsave( paste0(types[i],"survival_plot",".png"), plot = print(gg))
  #ggsave(filename=paste0(types[i],"survival_plot",".png"), plot=gg)
  dev.off()
  return(gg)
}