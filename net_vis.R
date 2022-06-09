net_vis <- function(data)
  
{
  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )
  print("Plotting the network...")
  list.of.packages <- c("devtools","dplyr","ggplot2","ggrepel","ggpubr","viridis","tibble","stringr",
                        "corrplot","tidyverse","igraph","visNetwork", "data.table", "CARNIVAL",
                        "viper", "CellNOptR","edgeR", "OmnipathR", "stringi","openxlsx","samr", "sna", "gplots","ggfortify","limma")
  
  invisible(lapply(list.of.packages, library, character.only = TRUE))
  
  
  signatures <- names(data)
  miRNAs <- c()
  total_rows <- 0 
  # first create the dataframe for the graph edges, vertices :
  
  for (i in 1:length(names(data)))
  {
    temp <- row.names(as.data.frame(data[[i]]$Table1))
    total_rows <- total_rows + nrow(data[[i]]$Table1) + nrow(data[[i]]$Table2)
    
    miRNAs <- c(miRNAs,temp)
    
  }
  df <- matrix(data = , nrow = total_rows, ncol = 3) 
  
  colnames(df) <- c("Signature","miRNA","Sign")
  
  
  
  
  miRNAs <- unique(miRNAs)
  
  index <- 1
  for (i in 1:length(names(data)))
  {
    miRNA_list <- row.names(as.data.frame(data[[i]]$Table1))
    for (j in 1:nrow(data[[i]]$Table1))
    {
      
      df[index, 1] <- names(data)[i]
      df[index, 2] <- miRNA_list[j]
      df[index, 3] <- 1
      index <- index + 1
      
      
    }
    miRNA_list <- row.names(as.data.frame(data[[i]]$Table2))
    for (jj in 1:nrow(data[[i]]$Table2))
    {
      
      df[index, 1] <- names(data)[i]
      df[index, 2] <- miRNA_list[jj]
      df[index, 3] <- -1
      index <- index + 1
      
      
    }
  }
  
  # write.csv(df,"df.csv")
  df <- as.data.frame(df)
  
  g <- igraph::graph.data.frame(df[,c(grep("Signature", colnames(df)),
                                      grep("miRNA", colnames(df)))])
  bipartite.mapping(g)
  
  V(g)$type <- bipartite_mapping(g)$type
  
  V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
  V(g)$shape <- ifelse(V(g)$type, "circle", "square")
  
  
  E(g)$color <- "lightgray"
  
  V(g)$label <- V(g)$name
  V(g)$degree <- igraph::degree(g)
  V(g)$size <- igraph::degree(g)
  
  
  nodes <- data.frame(id = V(g)$name, shape = ifelse(V(g)$type, "circle", "database"),color = ifelse(V(g)$type,"grey","green"))
  
  #shape = c("square", "triangle", "box", "circle", "dot", "star",
  # "ellipse", "database", "text", "diamond"),
  #nodes <- nodes[order(nodes$id, decreasing = F),]
  edges <- data.frame(from = df$Signature,
                      to = df$miRNA, color = ifelse(df$Sign==1,"blue","red"),
                      label = df$Sign, title = as.character(df[1:nrow(df),c(grep("Sign",
                                                                                 colnames(df)))]))
  
  
  
  
  mainP = paste0("From signatures to miRNAs and their positive/negative association")
  
  sp_g_5 <- visNetwork(nodes, edges, main = mainP,height = "1080px", width = "1920px") %>% visEdges(labelHighlightBold= "TRUE") %>%
    visInteraction(zoomView = TRUE) %>%
    visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)  %>% visLayout(improvedLayout = TRUE)
  
  
  
  visSave(sp_g_5, file = "FMB.html", background = "white")
  
  
  print("done!")
  
  
}