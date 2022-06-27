#CARNIVAL module

Carnival_opt <-function(iterator_index,df_EM,
                        results_dir,
                        inputs_dir,
                        disease_filename_j,
                        violin_gene,
                        cloned,
                        GAP,
                        cpu_threads,
                        network_similarity_threshold,
                        top_user,
                        hotspots,
                        top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user,load_other,radar_plot_data) {
  
  library(progress)
  library(dorothea)
  library(progeny)
  cplex_path = "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex"
  #cplex_path = "C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/bin/x64_win64"
  #  PREPARE CARNIVAL INPUT FILES:
  print("Preparing the viper regulon list...")
  write.csv(cloned,"cloned.csv")
  #df_EM <- df_EM[,-1]
  print(paste0("Standarizing Expression Matrix for the regulons of ", violin_gene, " ..."))
  
  df_EM <- scale(df_EM)
  df_EM <- as.data.frame(df_EM)
  
  # CREATE THE INPUT REGULON FOR CARNIVAL
  print("Creating the VIPER regulon object for CARNIVAL...")
  
  #regulon_df<- import_Omnipath_Interactions(select_organism = 9606)
  #regulon_df <- import_TFregulons_Interactions(select_organism = 9606) #keeping human regulons
  #save(regulon_df, file = "omni.RData")
  #=======================================================================
  # download all interactions from omnipath to build our full PNN
  
  # regulon_df <- import_Omnipath_Interactions(filter_databases=c())
  # regulon_df<- import_Omnipath_Interactions(select_organism = 9606)
  # temp_dir = getwd()
  # setwd(inputs_dir)
  # save(regulon_df, file = "omni.RData")
  # setwd(temp_dir)
  #=======================================================================
  load(file = paste0(inputs_dir,"/omni.RData"))
  
  #regulon_df <- regulons
  #regulon_df <- regulon_df[which(regulon_df$tfregulons_level%in%c("A", "B", "C", "D", "E")), ]
  #regulon_df <- regulon_df[which(regulon_df$tfregulons_level%in%c("A", "B", "C")), ]
  
  
  regulon_df <- regulon_df[which(regulon_df$is_directed==1), ] #keeping only directed
  regulon_df <- regulon_df[which((regulon_df$is_stimulation+regulon_df$is_inhibition)==1), ] #keeping only regulons which are either activations/inhibitions
  
  
  #regulon_df <- regulon_df[which(regulon_df$source_genesymbol %in% 
  #                                 regulons_violin_gene+regulon_df$target_genesymbol%in% regulons_violin_gene),]
  
  
  #### CRITICAL ### #### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ###
  #regulon_df <- regulon_df  %>%  dplyr::filter(regulon_df$target_genesymbol%in% regulons_violin_gene)
  #### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####### CRITICAL ####
  
  df <- matrix(data = , nrow = nrow(regulon_df), ncol = 3) #creating the regulon dataframe for the createRegulonList function
  
  df[, 1] = regulon_df$source_genesymbol
  df[, 3] = regulon_df$target_genesymbol
  
  df[which(regulon_df$is_stimulation==1), 2] <- "1" #assigning correct sign
  df[which(regulon_df$is_inhibition==1), 2] <- "-1"
  
  
  colnames(df) = c("Source", "Sign", "Target")
  df <- as.data.frame(unique(df))
  df$Source = as.character(df$Source)
  df$Sign = as.numeric(as.character(df$Sign))
  df$Target = as.character(df$Target)
  
  network <- df
  write.csv(network,paste0(results_dir,"/interactions.csv"))
  #save(network, file="network.RData")
  
  # regulon_object = createRegulonList(regulon_table = df)
  # 
  # save(regulon_object, file = "regulon_object.RData")
  setwd(results_dir)
  #load(file = paste0(inputs_dir,"/regulon_A_B_C_D_E.RData"))
  load(file = paste0(inputs_dir,"/regulon_A_B_C.RData"))
  
  #regulon_A_B_C <- regulon_A_B_C[which(regulon_A_B_C$tfregulons_level%in%c("A")), ]
  #############################################################################
  dir.create(path = paste0(results_dir,"/measurements"))
  
  if (load_other == TRUE) {
    carnival_path_1 <- paste0(results_dir, "/",disease_filename_j)
  }
  else {
    carnival_path_1 <- paste0(results_dir)
  }
  print(carnival_path_1)
  #dir.create(path = carnival_path_1)                        
  carnival_path <- paste0(carnival_path_1,"/opt")
  dir.create(path = carnival_path)
  
  
  setwd(paste0(results_dir,"/measurements"))
  save(network, file="network.RData")
  tp_input_i <- -1
  tp_input_a <- 1
  input_df_i  <- as.data.frame(tp_input_i)
  input_df_a  <- as.data.frame(tp_input_a)
  colnames(input_df_i) <- violin_gene
  colnames(input_df_a) <- violin_gene
  save(input_df_i, file="inputs_i.RData")
  save(input_df_a, file="inputs_a.RData")
  
  TF_activities = as.data.frame(viper::viper(eset = df_EM, 
                                             regulon = regulon_A_B_C, nes = T, 
                                             method = 'none', minsize = 4, 
                                             eset.filter = F)) ##estimating tf activities with viper
  
  
  print("Saving the regulon activities from Viper/DoRothEA...")
  save(TF_activities, file="measurements.RData")
  
  tfList <- generateDataframeTF(df = TF_activities, top = top_user, 1:ncol(df_EM)) ##generating the list of measurements. selecting the top50 measurements.
  
  mm <- get("model_human_full", envir = .GlobalEnv)
  pathway_activities <- progeny(expr = as.matrix(x = TF_activities), scale = TRUE, top = 100, perm = 1000, organism = "Human")
  
  
  load(file = paste0(inputs_dir,"//progenyMembers.RData")) # can also be modified by the user
  
  weightObj <- assignPROGENyScores(progeny = pathway_activities, progenyMembers = progenyMembers, access_idx = 1:nrow(pathway_activities))
  
  
  ##for each sample we generate a directory which will contain the dot figure and the list of result from carnival
  
  ##running carnival for each sample and saving each of the results in a list
  resList <- list()
  list_of_graphs <- c() # the list to contain all the graphs generated for comparison reasons later on
  
  
  # now create the dataframe which will contain the scores of comparing each graph newtwork with all other
  
  graph_heatmap <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) 
  
  # create similar binary size matrix that contains 1 if both row and column are same types of mutations
  tfList_names <- names(tfList)
  graph_heatmap_ID <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare same type of mutation
  graph_heatmap_DEL <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare deleterious or not
  graph_heatmap_hotspot <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare hotspot or not
  
  top_similar_1 <- c() 
  top_similar_2 <- c()
  top_bc_network <-  c()
  top_bc_vertex <- c()
  top_similar_score <- c()
  #length(tfList)
  total_Bar <- length(tfList)
  pb <- progress_bar$new(total = total_Bar)
  
  print("optimising networks with CARNIVAL MILP...")

  for (ii in 1: length(tfList)){
    pb$tick()
    Sys.sleep(1 / total_Bar)
    temp_dir <- paste0(carnival_path, "/", names(tfList)[ii])
    print(temp_dir)
    dir.create(path = temp_dir)
    
    
    final_input <- NULL
    final_input <- ifelse(toupper(toString(cloned$isDeleterious[ii])) == "TRUE", TRUE, FALSE)
    setwd(temp_dir)
    print(temp_dir)
    if(length(list.files(temp_dir)) ==0){
      if (!(is.null(final_input) )){
        if (final_input == TRUE) {
          print(paste0(violin_gene," knockdown..."))
          try(
            
            #invisible(readline(prompt="Press [enter] to continue"))
            res <- runCARNIVAL(inputObj = input_df_i, measObj = tfList[[ii]], netObj = network,weightObj = weightObj[[ii]],
                               solverPath = cplex_path , solver = "cplex",
                               dir_name = temp_dir, mipGAP = GAP, threads = cpu_threads)
            #poolrelGAP = 0, mipGAP = GAP, poolIntensity = 1, poolReplace = 2, limitPop = 10, poolCap = 10,
            # timelimit = 3600, threads = cpu_threads))
          )
        }
        else {
          print(paste0(violin_gene," stimulation..."))
          try(
            
            
            #invisible(readline(prompt="Press [enter] to continue"))
            res <- runCARNIVAL(inputObj = input_df_a, measObj = tfList[[ii]], netObj = network,weightObj = weightObj[[ii]],
                               solverPath = cplex_path , solver = "cplex",
                               dir_name = temp_dir, mipGAP = GAP, threads = cpu_threads)
            
            #poolrelGAP = 0, mipGAP = GAP, poolIntensity = 1, poolReplace = 2, limitPop = 10, poolCap = 10,
            #timelimit = 3600, threads = cpu_threads))
          )
        }
      }
      
      
    }
    else {print("Network has already been optimized ; folder not empty! Continuing...")}
  }
  
  all_networks <- list()
  pb <- progress_bar$new(total = total_Bar)
  
  print("Read all possible networks as graphs...")
  for (i in 1:  length(tfList)) {
    
    pb$tick()
    Sys.sleep(1 / total_Bar)
    #  Sys.sleep(0.1)
    # setTkProgressBar(pb, i, label=paste( round(i/total_Bar*100, 0),
    #                                     "% done"))
    
    net_base <- NULL
    try(temp_dir <- paste0(carnival_path, "/", names(tfList)[i]))
    try(base_file <- paste0(temp_dir,"/network_solution.dot"))
    
    try(dot_object <- sna::read.dot(base_file))
    try(net_base <- igraph::graph.adjacency(dot_object, mode = "directed"))
    all_networks[[i]] <- net_base
    
    # LEIDEN COMMUNITY DETECTION ALGORITHM
    #########################################
    #  V(net_base)$label <- V(net_base)$name
    # V(net_base)$label  <-sapply(strsplit(as.character(V(net_base)$label), split=' [', fixed=TRUE),function(x) (x[1]))
    # partition <- leiden(net_base,
    #                     resolution_parameter = 0.5, 
    #                     max_comm_size = 50, 
    #                     partition_type = "ModularityVertexPartition")
    # tp <- table(partition)
    # print(tp)
    # node.cols <- brewer.pal(max(c(10, partition)),"Pastel1")[partition]
    # leiden_object <- plot(net_base, vertex.color = node.cols)
    # print(leiden_object)
    ########################################
    
    
    # ================ LOUVAIN COMMUNITY DETECTION =======================
    g <- graph_from_literal()
    if (is.null(net_base) == FALSE) {
      g <- as.undirected(net_base)
      sp_g <- cluster_louvain(g)
      V(g)$community <- sp_g$membership
      net_data <- toVisNetworkData(g)
      V(g)$label <- V(g)$name
      V(g)$label <-sapply(strsplit(as.character(V(g)$label), split=' [', fixed=TRUE),function(x) (x[1]))
      
      
      
      gg <- g
      V(gg)$name <-sapply(strsplit(as.character(V(gg)$name), split=' [', fixed=TRUE),function(x) (x[1]))
      #########SPLIT ACROSS COMMUNITIES AND GET CENTRALITY FOR EACH SUBGRAPH AS SIGNATURE
      sub_graphs <- c()
      sub_objects <- c()
      unique_communities <- unique(V(gg)$community)
      temp_vertex <- c()
      for (ie in 1: length(unique_communities)){
        print(paste0("Community :", as.character(ie)))
        OV <- which(V(gg)$community == unique_communities[ie])
        g_subgraph_temp <- induced_subgraph(gg, OV)
        temp_vertex0 <- which.max(igraph::betweenness(g_subgraph_temp))
        print("Best centrality:")
        print(temp_vertex0)
        temp_vertex <- c(temp_vertex,names(temp_vertex0))
        #  readline(prompt="Press [enter] to continue")
      }
      signature <- paste(temp_vertex, collapse = '_') # signature
      print("Signature: ")
      print(signature)
      #readline(prompt="Press [enter] to continue")
      #########SPLIT ACROSS COMMUNITIES AND GET CENTRALITY FOR EACH SUBGRAPH AS SIGNATURE
      
      top_bc_network <-  c(top_bc_network, paste0(names(tfList)[i]))
      top_bc_vertex <- c(top_bc_vertex,signature)
      nodes <- data.frame(id =  V(g)$name,
                          group = V(g)$community,
                          label = V(g)$label,
                          shape = ifelse(V(g) %in% c(violin_gene), "star", "circle"))
      
      edges <- data.frame(net_data$edges,shadow = as.list(rep(TRUE, gsize(g))))
      mainP = paste0("Optimized network for:: ",names(tfList)[i],"::" , violin_gene)
      
      sp_g_5 <- visNetwork(nodes,edges, main = mainP, height = "700px", width = "100%") %>%
        visEdges(labelHighlightBold= "TRUE") %>%
        visNodes(size = 20) %>%
        visInteraction(zoomView = TRUE) %>%
        visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE,selectedBy = "group") %>%
        visPhysics(stabilization = TRUE)  %>% visIgraphLayout()
      visSave(sp_g_5, file = paste0(paste0(carnival_path, "/", names(tfList)[i]),"/Community_Detection",".html"), background = "white")
    }
    # ===============================================================
    
    
  }
  all_bc <- cbind(top_bc_network,top_bc_vertex)
  write.csv(all_bc,paste0(carnival_path,"/Betweenness_vertex_",disease_filename_j,".csv"))
  
  pb <- progress_bar$new(total = total_Bar)
  print("Starting optimized network comparisons...")
  for (i in 1:  length(tfList)) {
    pb$tick()
    Sys.sleep(1 / total_Bar)
    net_base <- all_networks[[i]]
    
    if (is.null(net_base) == FALSE) {
      
      for (jj in 1: length(tfList)) {
        score <- 0
        
        graph_heatmap[i,jj] <- 0
        graph_heatmap_ID[i,jj] <- 0
        graph_heatmap_DEL[i,jj] <- 0
        graph_heatmap_hotspot[i,jj]  <- 0
        
        
        if (i>jj) # only lower symmetric is enough
        {
          
          net_temp <- all_networks[[jj]]
          
          if (is.null(net_temp) == FALSE) {
            
            g_sim <- igraph::graph.intersection(net_base, net_temp, byname = "auto", keep.all.vertices = FALSE)
            
            if (gsize(net_base) > gsize(net_temp)) {
              score <- gsize(g_sim)/gsize(net_base)
              
            }
            else {
              
              score <- gsize(g_sim)/gsize(net_temp)
            }
            graph_heatmap[i,jj] <- score
            
            if (score >= 0.5) {
              top_similar_1 <-  c(top_similar_1, paste0(names(tfList)[i]))
              top_similar_2 <-  c(top_similar_2, paste0(names(tfList)[jj]))
              top_similar_score <- c(top_similar_score,score)
              
            }
            
            graph_heatmap_ID[i,jj] <- ifelse(cloned$Variant_Classification[i] == cloned$Variant_Classification[jj],1,0) # binary comparison
            #starts <- 15
            #ends <- 15
            #starts <- abs(cloned$Start_position[i] -cloned$Start_position[jj])
            #ends <- abs(cloned$End_position[i] -cloned$End_position[jj])
            #graph_heatmap_ID[i,jj] <- ifelse( starts <= 12 & ends <= 12,1,0) # binary comparison
           # graph_heatmap_DEL[i,jj] <- ifelse(cloned$isDeleterious[i] == cloned$isDeleterious[jj],1,0) # binary comparison / check if deleterious or not
            graph_heatmap_DEL[i,jj] <- ifelse(cloned$isDeleterious[i] == cloned$isDeleterious[jj],1,0)
            
            #graph_heatmap_hotspot[i,jj] <- ifelse(cloned$Protein_Change[i] %in% hotspots & cloned$Protein_Change[jj] %in% hotspots,1,0)
            graph_heatmap_hotspot[i,jj] <- ifelse(cloned$Protein_Change[i] %in% hotspots & cloned$Protein_Change[jj] %in% hotspots,1,0)
            
          }
          else{score <- 0}
        }
      }
    }
    else{score <- 0} # if base network iterating is already null
  }
  print("Finished network comparisons...")
  colnames(graph_heatmap) <- names(tfList)
  row.names(graph_heatmap) <- names(tfList)
  
  
  #print(graph_heatmap)
  setwd(results_dir)
  total_perc <- c() # similar
  total_perc_filtered <- c()  # similar and same type of mutation
  total_perc_filtered_del <- c() # similar and same flag for deleterious
  total_perc_filtered_del_S <- c() # similar and same type of mutation and deleterious flag
  total_perc_filtered_hotspot <- c() # similar and from same hotspot
  
  write.csv(graph_heatmap,paste0(carnival_path,"/Networks_Similarity_Scores_",disease_filename_j,".csv"))
  
  # now transform heatamp matrix to only account for similarity scores that come from same mutation types
  # this is done through filtering with graph_heatmap_ID as follows:
  graph_heatmap_filtered <- graph_heatmap*as.vector(graph_heatmap_ID) # this zeroes target based on zeroes in heatmap_D
  graph_heatmap_filtered_del <- graph_heatmap*as.vector(graph_heatmap_DEL) # this zeroes target based on zeroes in heatmap_D
  graph_heatmap_filtered_hotspot <- graph_heatmap*as.vector(graph_heatmap_hotspot)
  
  colnames(graph_heatmap_filtered) <- names(tfList)
  row.names(graph_heatmap_filtered) <- names(tfList)
  colnames(graph_heatmap_filtered_del) <- names(tfList)
  row.names(graph_heatmap_filtered_del) <- names(tfList)
  colnames(graph_heatmap_filtered_hotspot) <- names(tfList)
  row.names(graph_heatmap_filtered_hotspot) <- names(tfList)
  
  total_no_nets <- nrow(which( graph_heatmap >0, arr.ind=TRUE)) 
  total_no_nets_ID <- nrow(which( graph_heatmap_filtered >0, arr.ind=TRUE)) 
  total_n_nets_DEL <- nrow(which( graph_heatmap_filtered_del >0, arr.ind=TRUE)) 
  total_n_nets_hotspot <- nrow(which( graph_heatmap_filtered_hotspot >0, arr.ind=TRUE)) 
  
  
  for (i in 1:length(network_similarity_threshold)){
    #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
    temp_score <- sum(graph_heatmap >= network_similarity_threshold[i])
    temp_score <- round(temp_score/total_no_nets, digits = 2)
    total_perc <- c(total_perc,temp_score*100)
    
  }
  
  # now the same but only for same type of mutation:
  for (i in 1:length(network_similarity_threshold)){
    #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
    temp_score <- sum(graph_heatmap_filtered >= network_similarity_threshold[i])
    temp_score <- round(temp_score/total_no_nets_ID, digits = 2)
    total_perc_filtered <- c(total_perc_filtered,temp_score*100)
    
  }
  
  
  # now the same but only for same type of mutation:
  for (i in 1:length(network_similarity_threshold)){
    #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
    temp_score <- sum(graph_heatmap_filtered_del >= network_similarity_threshold[i])
    temp_score <- round(temp_score/total_n_nets_DEL, digits = 2)
    total_perc_filtered_del <- c(total_perc_filtered_del,temp_score*100)
    
  }
  
  # now the same but only for hotspots of TP53:
  for (i in 1:length(network_similarity_threshold)){
    #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
    temp_score <- sum(graph_heatmap_filtered_hotspot >= network_similarity_threshold[i])
    temp_score <- round(temp_score/total_n_nets_hotspot, digits = 2)
    total_perc_filtered_hotspot <- c(total_perc_filtered_hotspot,temp_score*100)
    
  }
  
  
  
  top_all <- cbind(top_similar_1,top_similar_2,top_similar_score)
  write.csv(top_all,paste0(carnival_path,"/top_similar_nets.csv"))
  #print(sp_CARNIVAL)
  write.csv(graph_heatmap_ID,paste0(carnival_path,"/graph_heatmap_ID.csv"))
  write.csv(graph_heatmap_DEL,paste0(carnival_path,"/graph_heatmap_DEL.csv"))
  write.csv(graph_heatmap_hotspot,paste0(carnival_path,"/graph_heatmap_hotspot.csv"))
  
  
  # general similarity
  # 
  # 
  
  ylimit <- length(tfList)
  xlimit <- ylimit
  if (xlimit > 100){ 
    xlimit <- 100
    ylimit <- xlimit
  }
  
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "%"," : ", 
                  paste(paste0(total_perc,"%"), collapse = ', '), " in ", disease_filename_j)
  sp_CARNIVAL <- (corrplot(graph_heatmap[1:xlimit,1:ylimit],type = "lower", col = cm.colors(100), title = title2,mar=c(0,0,1,0)))
  
  
  # same mutation similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and same type of mutation"," : ", 
                  paste(paste0(total_perc_filtered,"%"), collapse = ', '), " in ", disease_filename_j)
  sp_CARNIVAL <- (corrplot(graph_heatmap_filtered[1:xlimit,1:ylimit],type = "lower", col = cm.colors(100), title = title2,mar=c(0,0,1,0)))
  
  
  # deleterious or not similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and same deleterious flag"," : ", 
                  paste(paste0(total_perc_filtered_del,"%"), collapse = ', '), " in ", disease_filename_j)
  sp_CARNIVAL <- (corrplot(graph_heatmap_filtered_del[1:xlimit,1:ylimit],type = "lower", col = cm.colors(100), title = title2,mar=c(0,0,1,0)))
  
  
  # hotspot or not similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and same hotspot"," : ", 
                  paste(paste0(total_perc_filtered_hotspot,"%"), collapse = ', '), " in ", disease_filename_j)
  sp_CARNIVAL <- (corrplot(graph_heatmap_filtered_hotspot[1:xlimit,1:ylimit],type = "lower", col = cm.colors(100), title = title2,mar=c(0,0,1,0)))
  
    radar_plot_data[iterator_index,1] <- disease_filename_j
    radar_plot_data[iterator_index,2] <- total_perc[1]
    radar_plot_data[iterator_index,3] <- total_perc[2]
    radar_plot_data[iterator_index,4] <- total_perc[3]
    radar_plot_data[iterator_index,5] <- total_perc[4]
    
    radar_plot_data[iterator_index,6] <- total_perc_filtered[1]
    radar_plot_data[iterator_index,7] <- total_perc_filtered[2]
    radar_plot_data[iterator_index,8] <- total_perc_filtered[3]
    radar_plot_data[iterator_index,9] <- total_perc_filtered[4]
    
    radar_plot_data[iterator_index,10] <- total_perc_filtered_del[1]
    radar_plot_data[iterator_index,11] <- total_perc_filtered_del[2]
    radar_plot_data[iterator_index,12] <- total_perc_filtered_del[3]
    radar_plot_data[iterator_index,13] <- total_perc_filtered_del[4]
    
    radar_plot_data[iterator_index,14] <- total_perc_filtered_hotspot[1]
    radar_plot_data[iterator_index,15] <- total_perc_filtered_hotspot[2]
    radar_plot_data[iterator_index,16] <- total_perc_filtered_hotspot[3]
    radar_plot_data[iterator_index,17] <- total_perc_filtered_hotspot[4]
    
    print(radar_plot_data)
    write.csv(cloned,paste0(disease_filename_j,".csv"))
    write.csv(radar_plot_data,"radar_plot_data_temp.csv")
  return_list<- list("radar_plot_data" = radar_plot_data, sp_CARNIVAL,total_perc,total_perc_filtered,total_perc_filtered_del,total_perc_filtered_hotspot)
  
  return(return_list)
  
}