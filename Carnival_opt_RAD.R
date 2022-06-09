#CARNIVAL module

Carnival_opt_RAD <-function(df_EM,
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
                        top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user) {
  
  
  library(dorothea)
  library(progeny)
  
  #  PREPARE CARNIVAL INPUT FILES:
  print("Radiation data optimization...")
  #invisible(readline(prompt="Press [enter] to continue"))
  print("Preparing the viper regulon list...")
  
  df_EM <- df_EM[,-1]
  df_EM<- read.table(paste0(inputs_dir,"/radiation_fmb3.txt"),header=TRUE,sep="\t",row.names=1)
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
  load(file = paste0(inputs_dir,"/regulon_A_B_C_D_E.RData"))
  #load(file = paste0(inputs_dir,"/regulon_A_B_C.RData"))
  
  #regulon_A_B_C <- regulon_A_B_C[which(regulon_A_B_C$tfregulons_level%in%c("A")), ]
  #############################################################################
  dir.create(path = paste0(results_dir,"/measurements"))
  carnival_path_1 <- paste0(results_dir,"/opt")
  dir.create(path = carnival_path_1)                        
  carnival_path <- paste0(carnival_path_1,"/",disease_filename_j)
  dir.create(path = carnival_path)
  
  
  setwd(paste0(results_dir,"/measurements"))
  save(network, file="network.RData")
  tp_input_i <- NaN
  tp_input_a <- NaN
  input_df_i  <- as.data.frame(tp_input_i)
  input_df_a  <- as.data.frame(tp_input_a)
  colnames(input_df_i) <- violin_gene
  colnames(input_df_a) <- violin_gene
  save(input_df_i, file="inputs_i.RData")
  save(input_df_a, file="inputs_a.RData")
  
  #save(df_EM, file = "toy.RData")
  print(df_EM)
  TF_activities = as.data.frame(viper::viper(eset = df_EM, 
                                             regulon = regulon_A_B_C_D_E, nes = T, 
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
  top_similar_score <- c()
  
  total_Bar <- length(tfList)
  pb <- progress_bar$new(total = total_Bar)
 # pb <- winProgressBar(title = "progress bar", min = 0,max = length(tfList),width = 300)
  for (ii in 1: length(tfList)){
  #  setWinProgressBar(pb,ii,title = paste0(round(ii / length(tfList) * 100, 0),  paste0("% Optimizing networks in ",disease_filename_j)))
    pb$tick()
    Sys.sleep(1 / total_Bar)
    temp_dir <- paste0(carnival_path, "/", names(tfList)[ii])
    dir.create(path = temp_dir)
    
    final_input <- NULL
    
    # example of radiation condition : A2.R2_H460_2_PARENT_2_H460_PARENT
    
    # we split on hours to eithr activate or de-activate the perturbed gene:
    # 
    #final_input_s <- str_split(names(tfList)[ii], "_")

    final_input <- ifelse(str_detect(names(tfList)[ii],"H460_0_"), TRUE, FALSE)
    
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
  #close(pb)
  # length(tfList)
  # now compare the networks pair-wise
  genes_per_network <- c()
  
  #pb <- winProgressBar(title = "progress bar", min = 0,max = length(tfList),width = 300)
  for (i in 1:  length(tfList)) {
   # setWinProgressBar(pb,i,title = paste0(round(i / length(tfList) * 100, 0),  paste0("% Comparing optimized networks in ",disease_filename_j)))
    net_base <- NULL
    temp_dir <- paste0(carnival_path, "/", names(tfList)[i])
    base_file <- paste0(temp_dir,"/network_solution.dot")
    #g <- graph_from_literal()
    try(dot_object <- sna::read.dot(base_file))
    try(net_base <- igraph::graph.adjacency(dot_object, mode = "directed"))
    
    try(g <- net_base)
    
    #dot_lang_net <- readChar(base_file, file.info(base_file)$size)
    #visNetwork(dot = dot_lang_net,width = "100%")
    
    ###### statistics per network ########
    # print("###############################################")
    print(paste0("The optimized network contains: ", vcount(g), " nodes and ", ecount(g), " edges." ))
    vertex_list <- sapply(strsplit(as.character(names(V(g))), split=' [', fixed=TRUE),function(x) (x[1]))
    genes_per_network <- cbind(genes_per_network,vertex_list)
    
    # edge_list <- lapply(strsplit(as.character(E(g)), split=' [', fixed=TRUE),function(x) (x[1]))
    # 
    print(paste0("The network lacks ",round(length(setdiff(genes_HGNC_bkp,vertex_list))/length(genes_HGNC_bkp)*100,2) ,"% of MAPK genes."))
    print(paste0("The network lacks ",round(length(setdiff(regulons_violin_gene,vertex_list))/length(regulons_violin_gene)*100,2),"% of ",violin_gene," regulon."))
    # 
    # 
    # print("###############################################")
    #sp_g <- cluster_louvain(as.undirected(g))
    #V(g)$community <- sp_g$membership
    
    # net_data <- toVisNetworkData(g)
    # V(g)$label <- V(g)$name
    # V(g)$label <-sapply(strsplit(as.character(V(g)$label), split=' [', fixed=TRUE),function(x) (x[1]))
    #nodes <- data.frame(id =  V(g)$name, 
    #                     group = V(g)$community, 
    #                     label = V(g)$label,
    #                     color = ifelse("color=red" %in% as.character(V(g)$name),"red","yellow"),
    #                     shape = ifelse(V(g) %in% genes_HGNC_bkp, "star", "circle"))
    
     
    # edges <- data.frame(net_data$edges,shadow = as.list(rep(TRUE, gsize(g))), color = ifelse("color=red" %in% as.character(E(g)$name), "red", "blue"))
    # mainP = paste0("Optimized network for:: ",names(tfList)[i],":: mutation of " , violin_gene)
    # sp_g_5 <- visNetwork(nodes,edges, main = mainP, height = "1080px", width = "1920px") %>% visEdges(labelHighlightBold= "TRUE") %>%
    #   visInteraction(zoomView = TRUE) %>%
    #   visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
    #   visPhysics(stabilization = TRUE)  %>% visLayout(improvedLayout = TRUE)
    # visSave(sp_g_5, file = paste0(paste0(carnival_path, "/", names(tfList)[i]),"/Opt_Net_",disease_filename_j,".html"), background = "white")
    
    #invisible(readline(prompt="Press [enter] to continue"))
    
    if (is.null(net_base) == FALSE) {
      
      for (jj in 1: length(tfList)) {
        score <- 0
        net_temp <- NULL
        temp_dir <- paste0(carnival_path, "/", names(tfList)[jj])
        temp_file <- paste0(temp_dir,"/network_solution.dot")
        
        try(dot_object_temp <- sna::read.dot(temp_file))
        try(net_temp <- igraph::graph.adjacency(dot_object_temp))
        
        
        if (is.null(net_temp) == FALSE) {
          
          g_sim <- igraph::graph.intersection(net_base, net_temp, byname = "auto", keep.all.vertices = FALSE)
          
          if (gsize(net_base) > gsize(net_temp)) {
            score <- gsize(g_sim)/gsize(net_base)
            
          }
          else {
            
            score <- gsize(g_sim)/gsize(net_temp)
          }
          
        }
        
        # now store in the dataframe graph_heatmap the scores for net2net comparison results in similarity
        
        
        if (i<jj) # only lower symmetric is enough
        {
          graph_heatmap[i,jj] <- score
          
          if (score >= top_score) {
            top_similar_1 <-  c(top_similar_1, paste0(names(tfList)[i]))
            top_similar_2 <-  c(top_similar_2, paste0(names(tfList)[jj]))
            top_similar_score <- c(top_similar_score,score)
            
          }

          
         # graph_heatmap_ID[i,jj] <- ifelse(cloned$Variant_Classification[i] == cloned$Variant_Classification[jj],1,0) # binary comparison
          starts <- 15
          ends <- 15
          starts <- abs(cloned$Start_position[i] -cloned$Start_position[jj])
          ends <- abs(cloned$End_position[i] -cloned$End_position[jj])
          
          str_detect(names(tfList)[i],"H460_50B")
          str_detect(names(tfList)[jj],"H460_50B")
          
          graph_heatmap_ID[i,jj] <- ifelse(str_detect(names(tfList)[i],"H460_50B")
                                               && str_detect(names(tfList)[jj],"H460_50B") || str_detect(names(tfList)[i],"H460_60A")
                                                                                       && str_detect(names(tfList)[jj],"H460_60A"),1,0)


          graph_heatmap_DEL[i,jj] <- ifelse(str_detect(names(tfList)[i],"_H460_0_") == TRUE && str_detect(names(tfList)[jj],"_H460_0_") == FALSE,1,0)
          
          graph_heatmap_hotspot[i,jj] <- ifelse(str_detect(names(tfList)[i],"PARENT")
                                                && str_detect(names(tfList)[jj],"PARENT") || str_detect(names(tfList)[i],"RESISTANT")
                                                && str_detect(names(tfList)[jj],"RESISTANT"),1,0)
        }
        else 
        {
          graph_heatmap[i,jj] <- 0
          graph_heatmap_ID[i,jj] <- 0
          graph_heatmap_DEL[i,jj] <- 0
          graph_heatmap_hotspot[i,jj]  <- 0
        }
        
        
        
      }
    }
    else{score <- 0} # if base network iterating is already null
  }
  
  
 # close(pb)
  colnames(graph_heatmap) <- names(tfList)
  row.names(graph_heatmap) <- names(tfList)
  
  
  #print(graph_heatmap)
  setwd(results_dir)
  
  colnames(genes_per_network) <- names(tfList)
  write.csv(genes_per_network,"genes_per_net.csv")
  
  
  total_perc <- c() # similar
  total_perc_filtered <- c()  # similar and same type of mutation
  total_perc_filtered_del <- c() # similar and same flag for deleterious
  total_perc_filtered_del_S <- c() # similar and same type of mutation and deleterious flag
  total_perc_filtered_hotspot <- c() # similar and from same hotspot
  
  write.csv(graph_heatmap,paste0(carnival_path,"/Networks_Similarity_Scores_",disease_filename_j))
  
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
  write.csv(top_all,"top_similar_nets.csv")
  #print(sp_CARNIVAL)
  write.csv(graph_heatmap_ID,paste0(results_dir,"/graph_heatmap_ID.csv"))
  write.csv(graph_heatmap_DEL,paste0(results_dir,"/graph_heatmap_DEL.csv"))
  write.csv(graph_heatmap_hotspot,paste0(results_dir,"/graph_heatmap_hotspot.csv"))
  
  
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

  
  
  corr_matrix <- graph_heatmap
  colnames(corr_matrix) <- sapply(strsplit(colnames(corr_matrix), split='_', fixed=TRUE),function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  row.names(corr_matrix) <- sapply(strsplit(row.names(corr_matrix), split='_', fixed=TRUE), function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  
  
  sp_CARNIVAL <- (corrplot(corr_matrix,type = "upper", method="color", col = cm.colors(100), title = title2,mar=c(0,0,1,0),
                               tl.col = ifelse(str_detect(names(tfList),"PARENT"),"blue","red")))
  
  
  # same mutation similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and same cell line"," : ", 
                  paste(paste0(total_perc_filtered,"%"), collapse = ', '), " in ", disease_filename_j)
  
  
  corr_matrix <- graph_heatmap_filtered
  colnames(corr_matrix) <- sapply(strsplit(colnames(corr_matrix), split='_', fixed=TRUE),function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  row.names(corr_matrix) <- sapply(strsplit(row.names(corr_matrix), split='_', fixed=TRUE), function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  
  
  sp_CARNIVAL <- (corrplot(corr_matrix,type = "upper" , col = cm.colors(100), title = title2,mar=c(0,0,1,0),
                           tl.col = ifelse(str_detect(names(tfList),"PARENT"),"blue","red")))
  
  
  
  corr_matrix <- graph_heatmap_filtered_hotspot
  colnames(corr_matrix) <- sapply(strsplit(colnames(corr_matrix), split='_', fixed=TRUE),function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  row.names(corr_matrix) <- sapply(strsplit(row.names(corr_matrix), split='_', fixed=TRUE), function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  
  
  # deleterious or not similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and 0h to anything else"," : ", 
                  paste(paste0(total_perc_filtered_del,"%"), collapse = ', '), " in ", disease_filename_j)
  
  
  corr_matrix <- graph_heatmap_filtered_del
  colnames(corr_matrix) <- sapply(strsplit(colnames(corr_matrix), split='_', fixed=TRUE),function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  row.names(corr_matrix) <- sapply(strsplit(row.names(corr_matrix), split='_', fixed=TRUE), function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  
  
  sp_CARNIVAL <- (corrplot(corr_matrix,type = "upper" , col = cm.colors(100), title = title2,mar=c(0,0,1,0),
                           tl.col = ifelse(str_detect(names(tfList),"PARENT"),"blue","red")))
  
  

  # hotspot or not similarity
  title2 = paste0("Percentage of networks with similarity greater or equal than ",
                  paste(network_similarity_threshold*100, collapse = ', '), "% and same father"," : ", 
                  paste(paste0(total_perc_filtered_hotspot,"%"), collapse = ', '), " in ", disease_filename_j)
  
  
  corr_matrix <- graph_heatmap_filtered_hotspot
  colnames(corr_matrix) <- sapply(strsplit(colnames(corr_matrix), split='_', fixed=TRUE),function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  row.names(corr_matrix) <- sapply(strsplit(row.names(corr_matrix), split='_', fixed=TRUE), function(x) (paste(x[5],x[6],x[7],x[8],sep = "_")))
  
  
  sp_CARNIVAL <- (corrplot(corr_matrix,type = "upper" , col = cm.colors(100), title = title2,mar=c(0,0,1,0),
                           tl.col = ifelse(str_detect(names(tfList),"PARENT"),"blue","red")))
  
  return_list<- list(sp_CARNIVAL,total_perc,total_perc_filtered,total_perc_filtered_del,total_perc_filtered_hotspot)
  
  return(return_list)
  
}