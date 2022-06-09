library("gplots")
library("data.table")
mapk_data<-as.data.frame(read.csv("Carnival_EM.csv", row.names = 'Gene'), header = TRUE)
mapk_data <- mapk_data[,-1]
heat_clust <- as.data.frame(mapk_data)
row.names(heat_clust) <- mapk_data[,1]
heat_clust <- heat_clust[,-1]

# transpose

t_heat_clust <- data.table::transpose(heat_clust)

# get row and colnames in order
colnames(t_heat_clust) <- rownames(heat_clust)
rownames(t_heat_clust) <- colnames(heat_clust)


cc1 <- rainbow(nrow(t_heat_clust))
cc2 <- rainbow(ncol(t_heat_clust))



dend <-  heatmap.2(data.matrix(t_heat_clust), scale = "none", 
                   col = bluered(100), keysize = 1,
                   trace = "none", density.info = "none", 
                   key.title = "Log2 expression",
                   main = "Expression (log2) for the regulons VS mutation type",
                   xlab = paste0("Regulons of ", "violin_gene"), 
                   font.lab = 4,ylab = NULL, margins = c(6,24))