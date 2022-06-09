# dynamic graphs module


dynamic_graphs <-function(pairs_GC_CH_bkp,
                          violin_column,
                          violin_gene,
                          regulons_violin_gene,
                          dynamic_graphs_dir,
                          disease_filename_j,
                          skip_C){
###################################### IGRAPH PLOT ####################################
write.csv(pairs_GC_CH_bkp,"pairs_GC_CH_bkp.csv")
  
g <- graph.data.frame(pairs_GC_CH_bkp[,c(grep("Hugo_Symbol", colnames(pairs_GC_CH_bkp)),
                                         grep("Tumor_Sample_Barcode", colnames(pairs_GC_CH_bkp)))])
bipartite.mapping(g)

#####  OVERLAY VIOLIN CONNECTIONS FROM GRN to this cell-line->mutations graph:

V(g)$type <- bipartite_mapping(g)$type

V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
E(g)$color <- "lightgray"

targets <- NULL

if (!is.null(violin_column)) {
  targets <- incident(g, violin_gene)
}

V(g)$label <- V(g)$name
V(g)$degree <- igraph::degree(g)

if (!is.null(targets)) {
  E(g)$color[targets] <- 'green'
}
V(g)$size <- igraph::degree(g)

g_no_labels <- g

edge_attr(g) <- list(name = as.character(pairs_GC_CH_bkp[1:nrow(pairs_GC_CH_bkp),
                                                         c(grep("Variant_Classification", colnames(pairs_GC_CH_bkp)))]),
                     color = rep("red", gsize(g)))
edge_attr(g, "label") <- E(g)$name



main_G <- paste0("A bipartite graph of Genes->Cell Lines in ",
                 disease_filename_j, " \n Edges in green are ", violin_gene, " associated mutations only")
sp_g <- plot(g_no_labels, layout = layout.kamada.kawai, vertex.label.cex = 0.8,
             vertex.label.color = "black", edge.arrow.size=.5,main = main_G,
             rescale = FALSE, ylim=c(-5,7),xlim=c(-9,3), asp = 0 )




#SameLabel = function(e) {
#  V(g)[ends(g, e)[1]]$label == V(g)[ends(g, e)[2]]$label }
#g2 = delete_edges(g, which(sapply(E(g), SameLabel)))

########## modularity max ############
main_L <- paste0("Louvain community detection on a bipartite graph of Genes->Cell Lines in ", disease_filename_j)
sp_g_5 <- cluster_louvain(as.undirected(g))


V(g)$community <- sp_g_5$membership
nodes <- data.frame(id = V(g)$name, group = V(g)$community,
                    shape = ifelse(V(g)$type, "database", "circle"))

#shape = c("square", "triangle", "box", "circle", "dot", "star",
# "ellipse", "database", "text", "diamond"),
#nodes <- nodes[order(nodes$id, decreasing = F),]
edges <- data.frame(from = pairs_GC_CH_bkp[1:nrow(pairs_GC_CH_bkp),c(grep("Hugo_Symbol", colnames(pairs_GC_CH_bkp)))],
                    to = pairs_GC_CH_bkp[1:nrow(pairs_GC_CH_bkp),c(grep("Tumor_Sample_Barcode", colnames(pairs_GC_CH_bkp)))],
                    label = E(g)$name, title = as.character(pairs_GC_CH_bkp[1:nrow(pairs_GC_CH_bkp),c(grep("Protein_Change",
                                                                                                           colnames(pairs_GC_CH_bkp)))]),
                    shadow = as.list(rep(TRUE, gsize(g))))

mainP = paste0("A bipartite graph showing mutations from all genes (pathway and regulons) ",violin_gene," to cell lines")
sp_g_5 <- visNetwork(nodes, edges, main = mainP, height = "1080px", width = "1920px") %>% visEdges(labelHighlightBold= "TRUE") %>%
  
  visInteraction(zoomView = TRUE) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE)  %>% visLayout(improvedLayout = TRUE)


visSave(sp_g_5, file = paste0(dynamic_graphs_dir,"/Mutations_Network_",disease_filename_j,".html"), background = "white")

if (skip_C == FALSE) {
  # now save the same but only with the violin gene and its regulons
  TF_set <- c(regulons_violin_gene,violin_gene)
  TF_edges_df <- pairs_GC_CH_bkp %>% dplyr::filter(Hugo_Symbol %in% TF_set)
  g_TF <- graph.data.frame(TF_edges_df[,c(grep("Hugo_Symbol", colnames(TF_edges_df)),
                                          grep("Tumor_Sample_Barcode", colnames(TF_edges_df)))])
  bipartite.mapping(g_TF)
  V(g_TF)$type <- bipartite_mapping(g_TF)$type
  V(g_TF)$size <- igraph::degree(g)
  edge_attr(g_TF) <- list(name = as.character(TF_edges_df[1:nrow(TF_edges_df),
                                                          c(grep("Variant_Classification", colnames(TF_edges_df)))]),
                          color = rep("red", gsize(g_TF)))
  edge_attr(g_TF, "label") <- E(g_TF)$name
  
  sp_g_5 <- cluster_louvain(as.undirected(g_TF))
  
  V(g_TF)$community <- sp_g_5$membership
  
  TF_nodes <- data.frame(id = V(g_TF)$name, group = V(g_TF)$community,
                         shape = ifelse(V(g_TF)$type, "database", "circle"))
  
  
  TF_edges <- data.frame(from = TF_edges_df[1:nrow(TF_edges_df),c(grep("Hugo_Symbol", colnames(TF_edges_df)))],
                         to = TF_edges_df[1:nrow(TF_edges_df),c(grep("Tumor_Sample_Barcode", colnames(TF_edges_df)))],
                         label = E(g_TF)$name, title = as.character(TF_edges_df[1:nrow(TF_edges_df),c(grep("Protein_Change",
                                                                                                           colnames(TF_edges_df)))]),
                         shadow = as.list(rep(TRUE, gsize(g_TF))))
  
  mainTF = paste0("A bipartite graph showing mutations from the node of interest ",violin_gene," and its regulons to cell lines")
  sp_g_5_2 <- visNetwork(TF_nodes, TF_edges,main = mainTF,height = "1080px", width = "1920px") %>% visEdges(labelHighlightBold= "TRUE") %>%
    visInteraction(zoomView = TRUE) %>%
    visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = TRUE)  %>% visLayout(improvedLayout = TRUE)
  
  visSave(sp_g_5_2, file = paste0(dynamic_graphs_dir,"/Mutations_Network_TF_REGULONS_",disease_filename_j,".html"), background = "white")
  
}

return_list <- (sp_g)
return(return_list)
}