p53_final_plots <- function() {
  
  
  library(ggplot2)
  library(dplyr)

  p53_CCLE <- as.data.frame(read.csv("df_CCLE.csv"), header=TRUE)
  
  p53_TCGA <- as.data.frame(read.csv("df_TCGA.csv"), header=TRUE)
  
  
  
  
  p53_TCGA  <-  p53_TCGA %>% dplyr::group_by(Cancer_Type_full, status) 
  
  p53_TCGA <- p53_TCGA %>% dplyr::summarize(count = n())
  
  
  
  p53_CCLE  <-  p53_CCLE %>% dplyr::group_by(Site_Primary, status)
  
  p53_CCLE  <- p53_CCLE %>% dplyr::summarize(count = n())
  
  
  p53_CCLE <- p53_CCLE[!is.na(p53_CCLE$Site_Primary), ]
  
  p53_CCLE <- p53_CCLE %>% mutate(percent = count/sum(count))
  
  CCLE<- ggplot(p53_CCLE, aes(fill=status, y=count, x=Site_Primary)) +
    geom_bar(position="stack", stat="identity") +
    

    geom_text(aes(y = count, label = paste0(round(percent,digits = 2)*100, '%')),
              position = position_stack(vjust = 0.5), size = 4) +
    
    
    theme(axis.text.x=element_text(size=5, angle=90,hjust=0.95,vjust=0.2),plot.title = element_text(size = 9))+ font("xylab",size=10)+
    font("xy",size=10)+
    font("xy.text", size = 10) +
    font("legend.text",size = 10) + ggtitle("TP53 mutational frequency across cancers in CCLE")
  
  plot(CCLE) 
  
  #par(new=TRUE)
  
  
   
  p53_TCGA <- p53_TCGA[!is.na(p53_TCGA$Cancer_Type_full), ]
  #p53_TCGA <- add_column(p53_TCGA, Percentage = 0)
  p53_TCGA <- p53_TCGA %>% mutate(percent = count/sum(count))
  
  write.csv(p53_TCGA,"random.csv")
  
  
  TCGA<- ggplot(p53_TCGA, aes(fill=status, y=count, x=Cancer_Type_full)) + 
    geom_bar(position="stack", stat="identity") +
    

    geom_text(aes(y = count, label = paste0(round(percent,digits = 2)*100, '%')),
              position = position_stack(vjust = 0.5), size = 4) +
  
    theme(axis.text.x=element_text(size=5, angle=90,hjust=0.95,vjust=0.2),plot.title = element_text(size = 9))+ font("xylab",size=10)+
    font("xy",size=10)+
    font("xy.text", size = 10) +
    font("legend.text",size = 10) + ggtitle("TP53 mutational frequency across cancers in TCGA")
  
  plot(TCGA)

  

  return()
}