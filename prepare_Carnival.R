prepare_Carnival <- function(mapk_data,
                             filtered_expression_matrix,
                             results_dir,
                             inputs_dir,
                             disease_filename,
                             j,
                             violin_gene,
                             cloned,
                             GAP,
                             cpu_threads,
                             network_similarity_threshold, 
                             top_user,
                             hotspots,
                             top_score,
                             load_other,
                             TCGA_choice,
                             TCGA_disease,
                             do_GLM,
                             genes_e,
                             run_all_TCGA,
                             carnival_flag,
                             GLM_all,
                             rgenes_e,
                             eXML_dir,
                             reg_type,
                             key,genes_HGNC_bkp,regulons_violin_gene,FC_user,surv,FEM_user,ccle_iterator,condition)

{
  radar_plot_data <<- as.data.frame(matrix(data = , nrow = 1+length(TCGA_disease), ncol = 17))
  colnames(radar_plot_data) <- c("TCGA_SUBTYPE",	"Generic25",	"Generic50",	"Generic75",	"Generic90",	
                                 "MutationType25","MutationType50",	"MutationType75",	"MutationType90",	
                                 "Deleterious25",	"Deleterious50",	"Deleterious75",	"Deleterious90",	
                                 "Hotspot25",	"Hotspot50",	"Hotspot75",	"Hotspot90")
  
  # CCLE CASE:
  if (load_other == FALSE) {
    
    mapk_data_carnival <- merge(mapk_data,filtered_expression_matrix, by = "CELLLINE", all = TRUE)
    mapk_data_carnival$Variant_Classification[is.na(mapk_data_carnival$Variant_Classification)] <- "WT"
    setwd(results_dir)
    colnames(mapk_data_carnival) <- as.character(colnames(mapk_data_carnival))
    colnames(mapk_data_carnival) <- sapply(strsplit(colnames(mapk_data_carnival), split='..', fixed=TRUE),function(x) (x[1]))
    mapk_data_carnival <- mapk_data_carnival[,-(1:2),drop=F]
    mapk_data_carnival <- t(mapk_data_carnival)
    colnames(mapk_data_carnival) <- mapk_data_carnival[1,]
    mapk_data_carnival <- cbind(Gene = rownames(mapk_data_carnival), mapk_data_carnival)
    rownames(mapk_data_carnival) <- NULL
    mapk_data_carnival <- mapk_data_carnival[-1,,drop=F]
    mapk_data_carnival <- mapk_data_carnival[,colSums(is.na(mapk_data_carnival))<nrow(mapk_data_carnival)]
    
    
    
    
    write.csv(mapk_data_carnival,"Carnival_EM.csv")
    
    df_EM<-as.data.frame(read.csv("Carnival_EM.csv", row.names = 'Gene'), header = TRUE)
    
    
    if (FC_user == TRUE) {
      print("Using WT as control in CCLE...")
      # calculate Expression matrix using WT as control:
      df_EM_matrix <- data.matrix(df_EM)
      control_avg <- rowMeans(df_EM_matrix[,grep("WT", colnames(df_EM))])
      df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
      df_EM <- df_EM - control_avg
      
      
    }
    else {
      
      if (length(grep("WT", colnames(df_EM))) != 0) {
        df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
      }
      
    }
    
    write.csv(mapk_data_carnival,"Carnival_EM.csv")
    df_EM <- df_EM[,-1]
    
    if (condition=="Normal") {
      
      Carnival_opt_res <- Carnival_opt_RAD(df_EM,
                                           results_dir,
                                           inputs_dir,
                                           disease_filename[j],
                                           violin_gene,
                                           cloned,
                                           GAP,
                                           cpu_threads,
                                           network_similarity_threshold,
                                           top_user,
                                           hotspots,
                                           top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user)
      
    }
    else if (condition == "Radiation" ) {
      
      Carnival_opt_res <- Carnival_opt(ccle_iterator,df_EM,
                                       results_dir,
                                       inputs_dir,
                                       disease_filename[j],
                                       violin_gene,
                                       cloned,
                                       GAP,
                                       cpu_threads,
                                       network_similarity_threshold, 
                                       top_user,
                                       hotspots,
                                       top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user,load_other,radar_plot_data)
    }
    else if (condition == "Hypoxia_Large" ) {
      Carnival_opt_res <- Carnival_opt_HYPOXIA_L(df_EM,
                                                 results_dir,
                                                 inputs_dir,
                                                 disease_filename[j],
                                                 violin_gene,
                                                 cloned,
                                                 GAP,
                                                 cpu_threads,
                                                 network_similarity_threshold,
                                                 top_user,
                                                 hotspots,
                                                 top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user)
      
    }
    
    else if (condition == "Hypoxia_Small" ) {
      Carnival_opt_res <- Carnival_opt_HYPOXIA_S(df_EM,
                                                 results_dir,
                                                 inputs_dir,
                                                 disease_filename[j],
                                                 violin_gene,
                                                 cloned,
                                                 GAP,
                                                 cpu_threads,
                                                 network_similarity_threshold,
                                                 top_user,
                                                 hotspots,
                                                 top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user)
      
      
    }
    # CASE TCGA :
  }
  else 
  {
    
    for (i in 1:length(TCGA_disease)) {
      print("Starting iterator...")
      carnival_path_1 <- paste0(results_dir, "/", TCGA_disease[i])
      dir.create(path = carnival_path_1)  
      setwd(carnival_path_1)
      
      print(paste0("    Processing mutation and expresion data from TCGA ", TCGA_disease[i]," ..."))
      TCGA_TP53_mutations <- as.data.frame(read.table(file = paste0(inputs_dir,"/TCGA_mutations.tsv"), sep = '\t', header = TRUE))
      TCGA_TP53_mutations$Sample_ID <- sapply(substr(TCGA_TP53_mutations$Sample_ID, start = 1, stop = 12),function(x) (x[1]))
      TCGA_types <- as.data.frame(read.table(file = paste0(inputs_dir,"/subtype.txt"), sep = '\t', header = TRUE))
      isDeleterious <-  TCGA_TP53_mutations %>% dplyr::select(c("Sample_ID","Functional_Impact"))
      isDeleterious  <- add_column(isDeleterious , isDeleterious  = "False")
      isDeleterious <- isDeleterious %>% dplyr::mutate(isDeleterious = ifelse(stringr::str_detect(Functional_Impact, "impact: deleterious, score: 0"),"False","True"))
      colnames(isDeleterious) <- c("Sample_ID","Functional_Impact","isDeleterious")
      isDeleterious <-  isDeleterious %>% dplyr::select(c("Sample_ID","isDeleterious"))
      TCGA_TP53_mutations <- add_column(TCGA_TP53_mutations,isDeleterious = isDeleterious$isDeleterious)
      
      
      flag_out  = FALSE
      if (TCGA_choice == 1 || TCGA_choice == 3) {
        TCGA_disease_vector <- TCGA_types %>% dplyr::filter(Subtype %in% TCGA_disease[i])   
        
        TCGA_TP53_mutations <- TCGA_TP53_mutations %>% dplyr::filter(Sample_ID %in% TCGA_disease_vector$Patient_ID) 
        
      }
      else if (TCGA_choice == 2) {
        TCGA_TP53_mutations <- TCGA_TP53_mutations
      }
      if (nrow(TCGA_TP53_mutations)<2) {
        flag_out = TRUE
        next
      }
      write.csv(TCGA_TP53_mutations,"TCGA_TP53_mutations.csv")
      
      
      colnames(TCGA_TP53_mutations)[2] <- "CELLLINE"
      mapk_TCGA__data <-  TCGA_TP53_mutations %>% dplyr::select(c("CELLLINE","Mutation_Type","HGVSg","isDeleterious","Protein_Change"))
      colnames(mapk_TCGA__data) <- c("CELLLINE","Variant_Classification","Codon_Change","isDeleterious","Protein_Change")
      mapk_TCGA__data <- add_column(mapk_TCGA__data, Unique_Mut_ID = 0)
      mapk_TCGA__data <- mutate(mapk_TCGA__data, Unique_Mut_ID = sapply(mapk_TCGA__data$Variant_Classification,
                                                                        function(x) grep(x,unique(mapk_TCGA__data$Variant_Classification))))
      colnames(mapk_TCGA__data) <- c("CELLLINE", "Variant_Classification","Codon_Change","isDeleterious","Mutation_ID","Protein_Change")
      mapk_TCGA__data <- add_column(mapk_TCGA__data, Sample = mapk_TCGA__data$CELLLINE)
      non_unite_data <-  TCGA_TP53_mutations %>% dplyr::select(c("CELLLINE","Mutation_Type","HGVSg","isDeleterious","Protein_Change"))
      colnames(non_unite_data) <- c("CELLLINE","Variant_Classification","Codon_Change","isDeleterious","Protein_Change")
      mapk_TCGA__data <-mapk_TCGA__data %>% unite(Variant_Classification,Mutation_ID,Sample,Variant_Classification, Codon_Change, isDeleterious,Protein_Change, sep = "_", remove = TRUE)
      colnames(mapk_TCGA__data) <- c("CELLLINE", "Variant_Classification")
      
      
      load (paste0(inputs_dir,"/TCGA_expr.RData")) # gives back to workspace the TCGA expression matrix as TCGA_expr
      
      
      # now create the expression matrix for CARNIVAL from TCGA data:
      if (do_GLM == TRUE && GLM_all == TRUE && carnival_flag == FALSE) {
        
        filtered_expression_matrix_TCGA <- TCGA_expr
        
      }
      else {
        if (FEM_user == FALSE) {
          filtered_expression_matrix_TCGA <- TCGA_expr %>%  dplyr::select(matches(genes_e))
        }
        else {
          filtered_expression_matrix_TCGA <- TCGA_expr
        }
      }
      
      
      filtered_expression_matrix_TCGA$CELLLINE <- gsub("\\.","-",filtered_expression_matrix_TCGA$CELLLINE)
      filtered_expression_matrix_TCGA$CELLLINE <- sapply(substr(filtered_expression_matrix_TCGA$CELLLINE, start = 1, stop = 12),function(x) (x[1]))
      
      
      if (run_all_TCGA == TRUE) {
        filtered_expression_matrix_TCGA <-  filtered_expression_matrix_TCGA
      }
      else{
        filtered_expression_matrix_TCGA <-  filtered_expression_matrix_TCGA %>% dplyr::filter(CELLLINE %in% TCGA_disease_vector$Patient_ID) 
      }
      
      cloned <- merge(non_unite_data,filtered_expression_matrix_TCGA, by = "CELLLINE")
      
      
      genes_v <- paste0(violin_gene,"\\s*?\\.{2}ENSG000",collapse="|")
      genes_e_v <- paste0("^(", paste0("CELLLINE|",genes_v), ")")
      
      violin_gene_TCGA_expression <- filtered_expression_matrix_TCGA %>%  dplyr::select(matches(genes_e_v))
      
      colnames(violin_gene_TCGA_expression) <- c("CELLLINE","Expression_log2")
      
      
      
      cloned_violin <- merge(non_unite_data,violin_gene_TCGA_expression, by = "CELLLINE", all = TRUE)
      #####################
      WT_violin_TCGA <- filtered_expression_matrix_TCGA %>% dplyr::filter(!(CELLLINE %in% mapk_TCGA__data$CELLLINE))
      WT_violin_TCGA <- WT_violin_TCGA %>%  dplyr::select(c("CELLLINE"))
      mapk_TCGA__data2 <- merge(mapk_TCGA__data,WT_violin_TCGA,by = "CELLLINE" ,all = TRUE)
      
      
      
      
      cloned_violin <- merge(cloned_violin,WT_violin_TCGA, by = "CELLLINE",all = TRUE)
      cloned_violin <- add_column(cloned_violin, State = "MT")
      cloned_violin$State[is.na(cloned_violin$Variant_Classification)] <- "WT"
      cloned_violin$Variant_Classification[is.na(cloned_violin$Variant_Classification)] <- "WT"
      
      cloned_violin  <- mutate(cloned_violin, isHotspot = ifelse(Protein_Change %in% hotspots, TRUE, FALSE))
      write.csv(cloned_violin,"cloned_violin2.csv")
      
      # read TP53 CNA for TCGA: 
      calls_TP53_TCGA <- as.data.frame(read.table(file = paste0(inputs_dir,"/TP53_TCGA_CNA.txt"), sep = '\t', header = TRUE))
      calls_TP53_TCGA <-  calls_TP53_TCGA  %>% dplyr::select(c(2,3))
      colnames(calls_TP53_TCGA) <- c("CELLLINE", "GYSTIC_CALLS")
      calls_TP53_TCGA$CELLLINE <- sapply(substr(calls_TP53_TCGA$CELLLINE, start = 1, stop = 12),function(x) (x[1]))
      cloned_violin <- merge(cloned_violin,calls_TP53_TCGA, by = "CELLLINE")
      
      write.csv(cloned_violin,"cloned_violin.csv")
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels and types for ", violin_gene, " in TCGA ", TCGA_disease[i],
                     " \n (source data-sets: TCGA pancancer atlas ")
      
      
      sptcga77 <- ggplot(data = cloned_violin,aes(x = Variant_Classification, y = log2(Expression_log2+1), fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        #geom_point(shape = ifelse(total$CN_S == "Amplification", 8, 21),size=ifelse(total$CN_S == "Amplification", 3, 1),
        #           position = dodge, color= ifelse(total$CN_S == "Amplification", "red", "yellow"),alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in TCGA")) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) +
        stat_compare_means(method = "anova", label.y = 20, size = 8) +
        stat_compare_means(label.y = 24, size = 8)
      
      print(sptcga77)
      ggsave(filename="TCGA.png", plot=sptcga77)
      
      
      model <- aov(log2(Expression_log2+1)~Variant_Classification, data=cloned_violin)
      TM <- TukeyHSD(model, conf.level=.95)
      #GH <-games_howell_test(cloned_violin, log2(Expression_log2+1)~Variant_Classification, conf.level = 0.95, detailed = FALSE)
      
      print(TM)
      #print(GH)

      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in TCGA ", TCGA_disease[i])
      cloned_violin  <- mutate(cloned_violin, isHotspot = ifelse(Protein_Change %in% hotspots, TRUE, FALSE))
      write.csv(cloned_violin,"cloned_violin.csv")
      #cloned_violin[is.na(cloned_violin$Variant_Classification),] <- "WT"
      #cloned_violin <- cloned_violin[!is.na(cloned_violin$CN),]
      sptcga777 <- ggplot(data = cloned_violin,aes(x = State, y = log2(Expression_log2+1), fill = isHotspot))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in TCGA")) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) + 
        stat_compare_means(method = "anova", label.y = 20, size = 8)  + 
        stat_compare_means(label.y = 21, size = 8) +
        stat_compare_means(method = "t.test",label.y = 22, size = 8) +
        stat_compare_means(method = "wilcox.test",label.y = 23, size = 8)
      
      print(sptcga777)
      
      ggsave(filename="TCGA_WT_MT_isHotspot.png", plot=sptcga777)
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in TCGA ", TCGA_disease[i])
      cloned_violin  <- mutate(cloned_violin, isHotspot = ifelse(Protein_Change %in% hotspots, TRUE, FALSE))
      write.csv(cloned_violin,"cloned_violin.csv")
      #cloned_violin[is.na(cloned_violin$Variant_Classification),] <- "WT"
      #cloned_violin <- cloned_violin[!is.na(cloned_violin$CN),]
      sptcga777 <- ggplot(data = cloned_violin,aes(x = State, y = log2(Expression_log2+1), fill = isDeleterious))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in TCGA")) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) + 
        stat_compare_means(method = "anova", label.y = 20, size = 8)  + 
        stat_compare_means(label.y = 21, size = 8) +
        stat_compare_means(method = "t.test",label.y = 22, size = 8) +
        stat_compare_means(method = "wilcox.test",label.y = 23, size = 8)
      
      print(sptcga777)
      ggsave(filename="TCGA_WT_MT_isDel.png", plot=sptcga777)
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in TCGA ", TCGA_disease[i])
      cloned_violin  <- mutate(cloned_violin, isHotspot = ifelse(Protein_Change %in% hotspots, TRUE, FALSE))
      #write.csv(cloned_violin,"cloned_violin.csv")
      #cloned_violin[is.na(cloned_violin$Variant_Classification),] <- "WT"
      #cloned_violin <- cloned_violin[!is.na(cloned_violin$CN),]
      sptcga777 <- ggplot(data = cloned_violin,aes(x = State, y = log2(Expression_log2+1), fill = GYSTIC_CALLS))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = TRUE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in TCGA")) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) + 
        stat_compare_means(method = "anova", label.y = 20, size = 8)  + 
        stat_compare_means(label.y = 21, size = 8)
      
      print(sptcga777)
      ggsave(filename="TCGA_WT_MT_gystic.png", plot=sptcga777)
      
      ####################
      
      mapk_data_carnival <- merge(mapk_TCGA__data2,filtered_expression_matrix_TCGA, by = "CELLLINE")
      
      
      #####
      mapk_data_carnival <- mapk_data_carnival[order(mapk_data_carnival[,'Variant_Classification']), ]
      mapk_data_carnival$Variant_Classification[is.na(mapk_data_carnival$Variant_Classification)] <- "WT"
      mapk_data_carnival <- mutate(mapk_data_carnival,Variant_Classification = as.character(Variant_Classification))
      mapk_data_carnival[is.na(mapk_data_carnival)] <- 0
      #####
      #####
      surv_bkp <- mapk_data_carnival
      
      setwd(results_dir)
      
      colnames(mapk_data_carnival) <- as.character(colnames(mapk_data_carnival))
      colnames(mapk_data_carnival) <- sapply(strsplit(colnames(mapk_data_carnival), split='..', fixed=TRUE),function(x) (x[1]))
      mapk_data_carnival <- mapk_data_carnival[,-1,drop=F]
      #write.csv(mapk_data_carnival,"Carnival_EM_TCGA.csv")
      # if TCGA expression data was selected to do the GLM:
      
      if (do_GLM == TRUE) {
        
        print("    Using TCGA expression data for GLM...")
        filtered_expression_matrix_TCGA <- TCGA_expr
        filtered_expression_matrix_TCGA <- TCGA_expr %>%  dplyr::select(matches(rgenes_e))
        filtered_expression_matrix_TCGA$CELLLINE <- gsub("\\.","-",filtered_expression_matrix_TCGA$CELLLINE)
        filtered_expression_matrix_TCGA$CELLLINE <- sapply(substr(filtered_expression_matrix_TCGA$CELLLINE, start = 1, stop = 12),function(x) (x[1]))
        
        
        WT_violin_TCGA <- filtered_expression_matrix_TCGA %>% dplyr::filter(!(CELLLINE %in% mapk_TCGA__data$CELLLINE))
        WT_violin_TCGA <- WT_violin_TCGA %>%  dplyr::select(c("CELLLINE"))
        TCGA_GLM <- merge(mapk_TCGA__data,WT_violin_TCGA,by = "CELLLINE" ,all = TRUE)
        TCGA_GLM <- merge(TCGA_GLM,filtered_expression_matrix_TCGA, by = "CELLLINE")
        TCGA_GLM <- TCGA_GLM[order(TCGA_GLM[,'Variant_Classification']), ]
        TCGA_GLM$Variant_Classification[is.na(TCGA_GLM$Variant_Classification)] <- "WT"
        TCGA_GLM <- mutate(TCGA_GLM, Variant_Classification = as.character(Variant_Classification))
        TCGA_GLM[is.na(TCGA_GLM)] <- 0
        TCGA_GLM <- TCGA_GLM[,-1]
        
        
        if (key != "WT") {
          TCGA_GLM <- TCGA_GLM[-(grep("WT", TCGA_GLM$Variant_Classification)),]
        }
        
        saveRDS( TCGA_GLM, file = "regression.rds")
        write.csv(TCGA_GLM,"Renoir_TCGA.csv")
        
        resDir <- eXML_dir
        print("    Trying Generalized Linear Regression Modelling with TCGA ...")
        
        try(glregression <- MNR("regression.rds", 50, "TCGA",
                                reg_type, resDir, "alternative",key))
        stop()
        
      }
      if (carnival_flag == TRUE) {
        
        
        if (run_all_TCGA == TRUE) {
          TCGA_disease = c("TCGA_pancancer")
        }
        
        mapk_data_carnival <- t(mapk_data_carnival)
        colnames(mapk_data_carnival) <- mapk_data_carnival[1,]
        mapk_data_carnival <- cbind(Gene = rownames(mapk_data_carnival), mapk_data_carnival)
        rownames(mapk_data_carnival) <- NULL
        mapk_data_carnival <- mapk_data_carnival[-1,,drop=F]
        # remove full NA columns
        mapk_data_carnival <- mapk_data_carnival[,colSums(is.na(mapk_data_carnival))<nrow(mapk_data_carnival)]
        
        
        ##################################################################################
        print("    Using TCGA data for CARNIVAL...")
        
        write.csv(mapk_data_carnival,paste0("Carnival_EM_TCGA_",TCGA_disease[i],".csv"))
        #write.csv(mapk_data_carnival,"Carnival_EM_TCGA.csv")
        
        df_EM<-as.data.frame(read.csv(paste0("Carnival_EM_TCGA_",TCGA_disease[i],".csv"), row.names = 'Gene'), header = TRUE)
        
        
        
        df_EM_bkp <- df_EM
        
        if (FC_user == TRUE) {
          # calculate Expression matrix using WT as control:
          print("Using WT as control in TCGA...")
          df_EM_matrix <- data.matrix(df_EM)
          control_avg <- rowMeans(df_EM_matrix[,grep("WT", colnames(df_EM))])
          df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
          df_EM <- df_EM - control_avg
          
          
        }
        else {
          df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
          if (ncol(df_EM) <2) {df_EM <- df_EM_bkp}
        }
        df_EM <- df_EM[,-1]
        
        write.csv(df_EM,paste0("Carnival_EM_","TCGA_disease[i].csv"))
        
        
        try(
          Carnival_opt_res <- Carnival_opt(i,df_EM,
                                           results_dir,
                                           inputs_dir,
                                           TCGA_disease[i],
                                           violin_gene,
                                           cloned,
                                           GAP,
                                           cpu_threads,
                                           network_similarity_threshold, 
                                           top_user,
                                           hotspots,
                                           top_score,genes_HGNC_bkp,regulons_violin_gene,FC_user,load_other,radar_plot_data)
          
        )
        radar_plot_data <- Carnival_opt_res$radar_plot_data
        
        #print(radar_plot_data)
        #invisible(readline(prompt="Press [enter] to continue"))
        
        # process results
        sim_length <- length(network_similarity_threshold)
        carnival_csv <- as.data.frame(matrix(nrow=length(disease_filename)*(length(Carnival_opt_res)-1), ncol= 1 + sim_length))
        
        
        #############################################
        # for (iii  in 1:(length(Carnival_opt_res)-1)*j)
        # {
        #   carnival_csv[iii*j,1] <- disease_filename[j]
        #   for (ii in 2:(sim_length + 1)) 
        #   {
        #     colnames(carnival_csv)[1] <- "disease"
        #     colnames(carnival_csv)[ii] <- paste0(network_similarity_threshold[ii-1]*100,"%")
        #     
        #     carnival_csv[iii*j,ii] <- paste0(Carnival_opt_res[[iii+1]][ii-1],"%")
        #     
        #   }
        # }
        #############################################
        
        
        #  write.csv(carnival_csv,paste0(results_dir,"//carnival_summary.csv"))
      }
    }
    
    if (flag_out == FALSE) {
      print("Processing survival TCGA curated data...")
      opt_surv_samples <- c()
      opt_surv_DEL_status <- c()
      
      opt_surv_samples <- surv_bkp$CELLLINE
      opt_surv_DEL_status <- surv_bkp$Variant_Classification
      
      ################## create the survivability dataframe now ########################
      # opt_surv_samples <-  sapply(strsplit(as.character(colnames(mapk_data_carnival)), "_"),function(x) (x[1]))
      
      
      # the following defines which feature we will stratify upon the survival plot
      # opt_surv_DEL_status <- sapply(strsplit(as.character(colnames(mapk_data_carnival)), "_"),function(x) (x[2]))
      # opt_surv_DEL_status <- as.character(colnames(surv_bkp))
      # opt_surv_DEL_status <- sapply(strsplit(as.character(opt_surv_DEL_status), "_"),function(x) (x[length(x)]))
      
      
      
      
      df <- matrix(data = , nrow = length(opt_surv_samples), ncol = 2) 
      
      df[, 1] = opt_surv_samples
      df[, 2] = opt_surv_DEL_status
      
      
      
      colnames(df) = c("bcr_patient_barcode", "DEL_function")
      df <- as.data.frame(df)
      
      
      surv_2 <- merge(df,surv,by = "bcr_patient_barcode" , all = FALSE)
      
      write.csv(surv_2,"opt_surv.csv")
      print("Survivability data merged with optimization matrix")
      
      #sap <- NULL # survival analysis plot
      
      
      #sap <- surv_analysis()
      #print(sap)
      write.csv(radar_plot_data,"radar_plot_data_final.csv")
    }
    
    
  } #end of TCGA case
  
  
  
} # end of function