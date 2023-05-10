CCLE2 <-function(disease_name, 
                 target_gene,
                 do_GLM,reg_type, 
                 carnival_flag, 
                 merge_pathways,
                 only_regulons,
                 GLM_all,
                 new_version,
                 dataset_version,
                 GAP,
                 cpu_threads,
                 sig_user,
                 load_other,
                 top_user,
                 TCGA_disease,
                 run_all_TCGA,
                 hotspots,
                 TCGA_choice,
                 key,
                 top_score,
                 load_other_GLM,
                 GLM_signature,
                 FC_user,
                 FEM_user,
                 GLM_predict_user,condition)  {
  
  
  # SHINY VERSION 
  ###########################################################################################################
  # This function analyzes RNAseq data from the CCLE database and uses a user selected Gene Regulatory Network
  # to collect a list of genes and check whether a subset of those were mutated in a specific disease (user)
  # by reporting a pairwise Student T-test of the means of the expression profiles of the mutant cell lines
  # versus the non-mutant (abused terminology as WT) and the log2 fold change of them.
  # Copyright : Charalampos Triantafyllidis, Francesca M. Buffa
  #             Department of Oncology, University of Oxford, UK, 2019-2020
  # Output:
  #       pdf file with all plots collated for the disease as selected by user
  # Charalampos P. Triantafyllidis, 2020
  ############################################################################################################
  
  ######################### MODIFY AS REQUIRED BELOW ############################
  ###############################################################################
  #    USERS MUST DEFINE THE WD IN THE PARENT FOLDER WHERE THIS FUNCTION IS SAVED, MANUALLY
  #    THEY MUST ALSO DEFINE THE PATH BELOW FOR CPLEX , REQUIRED TO RUN CARNIVAL LATER ON
  options(warn=-1)
  wd = getwd()
  
  #cplex_path = "C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/bin/x64_win64"
  ###############################################################################
  # load the multinomial regression package
  # this predicts the type of mutation of the violin_gene based on
  # the expression levels of its regulons across all CCCLE samples
  source(paste0(wd,'//analysis.R'))
  source(paste0(wd,'//analysis_predict.R'))
  source(paste0(wd,'//runDoRothEA.R'))
  source(paste0(wd,'//createRegulonList.R'))
  source(paste0(wd,'//generateDataframeTF.R'))
  # CCLE modules:
  source(paste0(wd,'//pca_single.R'))
  
  source(paste0(wd,'//pancancer_graphs.R'))
  source(paste0(wd,'//dynamic_graphs.R'))
  source(paste0(wd,'//Carnival_opt.R'))
  source(paste0(wd,'//Carnival_opt_HYPOXIA.R'))
  source(paste0(wd,'//Carnival_opt_RAD.R'))
  source(paste0(wd,'//Carnival_opt_HYPOXIA_L.R'))
  source(paste0(wd,'//Carnival_opt_HYPOXIA_S.R'))
  source(paste0(wd,'//createRegulons.R'))
  source(paste0(wd,'//createRegulonList.R'))
  source(paste0(wd,'//prepare_GLM_data.R'))
  source(paste0(wd,'//prepare_Carnival.R'))
  source(paste0(wd,"//assignPROGENyScores.R"))
  source(paste0(wd,"//surv_analysis.R"))
  ###############################################################################
  
  # hotspot mutations (protein change in CCLE terms) :
  # https://doi.org/10.1038/cdd.2017.180
  hotspots_bkp <- c("p.R175H","p.R248Q","p.R273H","p.R248W", "p.R273C", "p.R282W", "p.G245S")
  
  ################################# OPTIONAL ################################################################
  # Reading the supplied .csv files into matrices
  
  # ######## READING THE CSV INPUT FILES FROM CCCLE DATABASE#########
  # print("Reading expression matrix...")
  # expr_matrix <- as.data.frame(read.csv("C:/Users/Harry/Desktop/DepMap19Q3/expressions.csv", header = TRUE))
  # print("Reading mutation matrix...")
  # mut_matrix  <- as.data.frame(read.csv("C:/Users/Harry/Desktop/DepMap19Q3/mutations.csv", header = TRUE))
  # print("Reading gene list from specific network...")
  # network_genes <- as.data.frame(read.csv("C:/Users/Harry/Desktop/DepMap19Q3/network_genes.csv", header = TRUE))
  # print("Reading copy number data...")
  # cn_csv <- as.data.frame(read.csv("C:/Users/Harry/Desktop/DepMap19Q3/CCLE_gene_cn.csv", header = TRUE))
  # print("Import of required data-sets done.")
  # 
  ################################# OPTIONAL ################################################################
  
  
  inputs_dir = paste0(wd,"//inputs")
  if (load_other == TRUE) {
    results_dir = paste0(wd,"//TCGA")
  }
  else {
    results_dir = paste0(wd,"//CCLE")
  }
  print(results_dir)
  dir.create(path = results_dir)
  results_dir_bkp <- results_dir
  
  if (new_version==TRUE) {
    # case new CCLE release is out and you need to directly read the new files from the input folder:
    new_version_folder <- as.character(dataset_version)
    print(paste0("Preparing input dataset for new CCLE version ", as.character(dataset_version[1])))
    print("Reading copy number matrix...")
    cn_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CCLE_gene_cn.csv"), header=TRUE))
    print("Reading expression matrix...")
    expr_matrix_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CCLE_expression_full.csv"), header=TRUE))
    print("Reading mutation profiles matrix...")
    mut_matrix_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CCLE_mutations.csv"), header=TRUE))
    print("Reading RNAseq reads matrix...")
    RNAseq_matrix <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CCLE_RNAseq_reads.csv"), header=TRUE))
    print("Reading RPPA data matrix...")
    RPPA <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CCLE_RPPA.csv"), header=TRUE))
    
   try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="Annotation_Transcript")] <- "Tumor_Sample_Barcode"
   )
      
    
    save(expr_matrix_csv,
         mut_matrix_csv, 
         cn_csv,
         RNAseq_matrix, 
         gene_effect,
         common_essentials,
         gene_dep,
         RPPA,file=paste0(inputs_dir,"/CCLE_",
                          as.character(dataset_version),".RData"))
    
    
    print("Now reading Achilles CRISPR data...")
    # reading CRISPR data:
    gene_effect <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/Achilles_gene_effect.csv"), header=TRUE))
    common_essentials <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/Achilles_common_essentials.csv"), header=TRUE))
    gene_dep <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/CRISPR_gene_dependency.csv"), header=TRUE))
    
    
    save(gene_effect,
         common_essentials,
         gene_dep,
         file=paste0(inputs_dir,"/CRISPR_CCLE_" , as.character(dataset_version),".RData"))
    
    
    # reading TCGA data
    print("Now reading RNaseq TCGA matrix...")
    
    try(TCGA_RNA <- read.table(file = paste0(inputs_dir,'/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'), sep = '\t', header = TRUE))
    print("Saving new TCGA RNAseq matrix...")
    try(save(TCGA_RNA,file=paste0(inputs_dir,"/TCGA_RNA",".RData")))
    
  }
  
  
  GLM_sig <- c()
  if (GLM_signature == TRUE) {
    print("Reading Regression Signature ...")
    GLM_sig <- as.data.frame(read.csv(paste0(inputs_dir,"/Report.csv"), header=TRUE))
    GLM_sig <- sapply(strsplit(GLM_sig[2:nrow(GLM_sig),1], split='..', fixed=TRUE),function(x) (x[1]))
    GLM_sig <- sapply(strsplit(GLM_sig, split='.', fixed=TRUE),function(x) (x[1]))
    write.csv(GLM_sig,paste0(results_dir,"/GLM_sig.csv"))
  }
  
  
  
  
  
  setwd(inputs_dir)
  print("Loading primary data-sets (CCLE) ...")
  load(file=paste0("CCLE_",dataset_version,".RData"))
  print("Done...")
  
  print("Checking column names for newer version of CCLE mutation profiles...")
  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="DepMap_ID")] <- "Tumor_Sample_Barcode"
      )
  print("Reading the cell line names mapping and suppliers list...")
  names_dir = paste0(inputs_dir, "//", "cell_line_names.csv")
  cell_line_names <- as.data.frame(read.csv(names_dir, header = TRUE))
  names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
  colnames(names_mat) <- c("CELLLINE","CCLE_ID")
  
  
  
  plot_list <- c()
  
  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )
  
  print("Loading libraries required...")
  list.of.packages <- c("dplyr","ggplot2","ggrepel","ggpubr","viridis","tibble","stringr",
                        "corrplot","tidyverse","igraph","visNetwork", "data.table", "CARNIVAL",
                        "viper", "CellNOptR", "OmnipathR", "stringi","openxlsx",
                        "sna", "gplots","ggfortify","limma", "UpSetR","survival", "survminer","ggcorrplot")# "edgeR",
  
  invisible(lapply(list.of.packages, library, character.only = TRUE))
  
  # user specifies the name of the cancer type to analyze (can input also "all" to analyze all of the supplied cell line files)
  disease_filename <- c()
  
  # measuring run-time
  start <- Sys.time()
  
  # TP53, CDC20, CENPA, KIF2C, PLK1 etc.
  violin_gene <- target_gene
  
  
  # rename some columns to tag cell lines in expression, mutation and CN matrices :
  colnames(RNAseq_matrix)[1] <- "CELLLINE"
  colnames(expr_matrix_csv)[1] <- "CELLLINE"
  colnames(cn_csv)[1] <- "cell_lines"
  colnames(RPPA)[1] <- "disease"
  
  
  vg <- paste0(violin_gene,"..ENSG000", collapse="|")
  violin_gene_E_f <- paste0("^(", paste0("CELLLINE|",vg), ")")
  violin_gene_E <- expr_matrix_csv %>% dplyr::select(matches(violin_gene_E_f))
  cell_line_dir = paste0(inputs_dir, "//", "CELL_LINES")
  
  
  if (disease_name == "all") {
    disease_filename <- list.files(cell_line_dir)
  }
  else {
    disease_filename <- c(disease_filename,paste0(disease_name , ".csv"))
  }
  # csv file to save the results from CARNIVAL for the similarity comparison of the networks
  network_similarity_threshold = c(0.25,0.5,0.75,0.9)
  
  
  print("Creating pdf results file...")
  setwd(results_dir)
  # create graph visual dir
  dynamic_graphs_dir = paste0(results_dir,"/dynamic_graphs")
  dir.create(path = dynamic_graphs_dir)
  
  eXML_dir <- paste0(results_dir,"/GLM")
  if (do_GLM == TRUE) {dir.create(path = eXML_dir)} 
  
  
  all_cancers <- NULL
  all_cancers_all_genes <- NULL
  
  # read the regulons of the violin_gene
  # first get the path for the csv that contains the regulons
  
  # tp53 sig from CELL REPORTS   https://doi.org/10.1016/j.celrep.2019.07.001
  tp53_sig <- c("CDC20","CENPA","KIF2C","PLK1") # 
  
  #tp53_sig <-c ()
  pnas_sig <- c()
  #pnas_sig <- c("MYBL2", "TFF1", "BRRN1",
  #              "CHAD", "SCGB3A1", "DACH","CDCA8","LAF4",
  #              "NY-BR-1", "DACH", "MYBL2", "CACNG4", "CYBRD1","LRP2", "SCGB3A1", "TFF1",
  #              "STC2", "AGR2") # https://doi.org/10.1073/pnas.0506230102
  
  print(paste0("Reading the regulons of ", violin_gene," from DoRothEA (all confidence levels)"))
  regulons_dir = paste0(inputs_dir, "//", "regulons.csv")
  regulons_csv <- as.data.frame(read.csv(regulons_dir, header = TRUE))
  
  print("Reading TCGA survival data...")
  surv <- as.data.frame(read.csv(paste0(inputs_dir,"/curated_surv.csv"), header=TRUE))
  print("Done!")
  
  # print("Leaving out regulons that came from ARACNe-GTEx only...")
  # now filter only for violin gene as source TF to get the list of its regulons - leave GTEx source out
  #regulons <- regulons_csv  %>%  dplyr:: filter(source_genesymbol %in% violin_gene & !(stringr::str_detect(sources, "^ARACNe-GTEx$")))
  regulons <- regulons_csv  %>%  dplyr:: filter(source_genesymbol %in% violin_gene)
  regulons_violin_gene <- NULL
  
  regulons_violin_gene <- unique(as.character(regulons$target_genesymbol))
  setwd(results_dir)
  # create violin expression data folder
  violin_expr_path = paste0(results_dir,"/", violin_gene, "_expression_data")
  violin_mut_path = paste0(results_dir,"/", violin_gene, "_mutation_data")
  cell_line_IDs_names = paste0(results_dir,"/", "_cell_line_IDs_names_data")
  pca_objects = paste0(results_dir,"/", "pca_data")
  dir.create(path = violin_expr_path)
  dir.create(path = violin_mut_path)
  dir.create(path = cell_line_IDs_names)
  dir.create(path = pca_objects)
  write.csv(regulons_violin_gene,paste0(violin_gene,"_regulons.csv"))
  setwd(inputs_dir)
  
  print("Reading the pathway regulatory network...")
  # read the GRN:
  GRN_dir = paste0(inputs_dir, "//", "GRN.graphml")
  GRN <- read_graph(GRN_dir, format = "graphml")
  
  genes_HGNC  = V(GRN)$name
  genes_HGNC_bkp <- genes_HGNC
  
  # define whether user regulons/sig or not
  if (is.null(sig_user) == FALSE) {
    regulons_violin_gene <- unique(c(sig_user))
    # genes_HGNC <- c(regulons_violin_gene, "TP53")
    # genes_HGNC_bkp <- genes_HGNC
  }
  
  if (merge_pathways == TRUE) {
    print("Merging pathway genes with the regulons of the target gene...")
    genes_HGNC <- unique(c(genes_HGNC,violin_gene,regulons_violin_gene))
  }
  
  
  # ******************************************************
  if (GLM_signature == TRUE)
  {
    genes_HGNC <- unique(c(GLM_sig,"TP53"))
    regulons_violin_gene <- GLM_sig
  }
  
  if (only_regulons == TRUE){
    
    #genes_HGNC <- unique(c(regulons_violin_gene))
    genes_HGNC <- unique(c(regulons_violin_gene,violin_gene))
  }
  # print set diff
  
  test_bed <- list(regulon_Dorothea = as.list(unique(as.character(regulons$target_genesymbol))), GLM_sig = as.list(GLM_sig),
                   CELL_REPORTS_TP53_sig = tp53_sig, MAPK_pathway_genes = as.list(V(GRN)$name) )
  print(upset(fromList(test_bed), order.by = "freq"))
  
  ################################################### MAIN LOOP ##################################################
  for (j in 1: length(disease_filename)) {
    ccle_iterator <- j
    if (load_other == TRUE) {
      results_dir <- results_dir_bkp
      setwd(results_dir)
    }
    else {
    temp_dir <- paste0(results_dir_bkp ,"/cancers_",disease_filename[j])
    dir.create(path = temp_dir)
    results_dir <- temp_dir
    setwd(temp_dir)
    }

    pdffile = paste0("Plots_cancers_",disease_filename[j],"_",dataset_version,".pdf")
    pdf(file = pdffile, height = 10, width = 18, onefile = TRUE)
    ###### separating title page for each cancer type #####
    a = paste0("Analyses for ", disease_filename[j])
    
    
    plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
         xaxt='n', yaxt='n', xlab='', ylab='')
    text(1,4,a, pos=4)
    command <- paste0(cell_line_dir,"//",disease_filename[j])
    command <- shQuote(command)
    command <- paste0("read.csv(",command)
    command <- paste0(command,",header = TRUE)")
    command <-paste0("as.data.frame(",command)
    command <- paste0(command,")")
    
    if (disease_name == "all_cell_lines") {
      disease_cell_lines = expr_matrix_csv$CELLLINE
    }
    else {
      disease_csv <- eval(parse(text = command))
      disease_cell_lines = disease_csv[, 1]
    }
    #violin_column <- rep(0, length(disease_cell_lines))
    violin_column <- NULL
    
    ################################Processing names of genes#################
    exact_genes <- NULL
    exact_genes <- paste0("^",genes_HGNC, "$", collapse="|")
    
    genes <- paste0(genes_HGNC,"\\s*?\\.{2}ENSG000",collapse="|")
    #genes <- paste0(genes_HGNC,"..ENSG000", collapse="|")
    genes0 <- paste0("\\b",genes_HGNC,"\\s*?\\.{2}\\b", collapse="|")
    genes_e <- paste0("^(", paste0("CELLLINE|",genes), ")")
    genes_cn <- paste0("^(", paste0("cell_lines|",genes0), ")")
    tumor_samples <- paste0("^(", paste(disease_cell_lines, collapse="|"), ")")
    
    print(paste0("Filtering the expression profiles in ",disease_filename[j]))
    
    ############################################################################################################
    filtered_RPPA <- RPPA  %>%   dplyr::filter(stringr::str_detect(disease, "_BREAST"))
    
    filtered_expression_matrix <- expr_matrix_csv  %>%  dplyr::filter(CELLLINE %in% disease_cell_lines)
    filtered_expression_matrix_disease <- filtered_expression_matrix
    # create a csv with the cell line names and BROAD IDs only for the specific disease
    names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
    colnames(names_mat) <- c("CELLLINE","CCLE_ID")
    
    small_filtered_expression_matrix <- filtered_expression_matrix %>% dplyr::select(c("CELLLINE"))
    cell_line_map <- merge(names_mat,small_filtered_expression_matrix, by = "CELLLINE")
    write.csv(cell_line_map,paste0(cell_line_IDs_names,"/","_cell_line_mapping_",disease_filename[j]))
    colnames(names_mat) <- c("cell_lines","CCLE_ID")
    # prepare the calls for the GYSTIC data for violin_gene
    print(paste0("Reading GYSTIC cBioPortal calls for ", violin_gene,"..."))
    calls <- read.table(file = paste0(inputs_dir,"/calls.txt"), sep = '\t', header = TRUE)
    colnames(calls)[2] <- "CCLE_ID"
    calls <- merge(calls,names_mat, by = "CCLE_ID")
    
    
    if (ncol(violin_gene_E)!=2) {
      print(paste0("No expression data was found for ", violin_gene, ". Exiting..."))
      stop()
    }
    else{
      setwd(results_dir)
      
      # selected expression data only for the violin gene:
      violin_expression_data <- NULL
      names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
      colnames(names_mat) <- c("CELLLINE","CCLE_ID")
      violin_expression_data <- filtered_expression_matrix %>%  dplyr::select(matches(violin_gene_E_f))
      violin_expression_data <- merge(names_mat,violin_expression_data, by = "CELLLINE")
      write.csv(violin_expression_data,paste0(violin_expr_path,"/",violin_gene,"_expression_in_",disease_filename[j]))
    }
    
    filtered_expression_matrix <- filtered_expression_matrix %>%  dplyr::select(matches(genes_e))
    
    if (ncol(filtered_expression_matrix) <=2) {
      print(paste0("No expression data was found for your signature. Exiting..."))
      stop()
    }
    
    setwd(results_dir)
    
    ############################################################################################################
    print(paste0("Filtering the mutation profiles across only the GRN and Regulons in ",disease_filename[j]))
    filtered_mutation_matrix <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines &
                                                                   !(Variant_Classification %in% "Silent"))
    
    skip_C = FALSE
    # we now filter only for the genes that also belong to our GRN
    filtered_mutation_matrix <- filtered_mutation_matrix %>% dplyr::filter(stringr::str_detect(Hugo_Symbol, exact_genes))
    if (!(violin_gene %in% filtered_mutation_matrix$Hugo_Symbol)) {
      print(paste0("No mutation data found for ", violin_gene, "!"))
      skip_C = TRUE
    }
    
    
    if (skip_C == FALSE) {
      setwd(results_dir)
      
      # selected expression data only for the violin gene:
      violin_mutation_data <- NULL
      names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
      colnames(names_mat) <- c("Tumor_Sample_Barcode","CCLE_ID")
      violin_mutation_data <- filtered_mutation_matrix %>%  dplyr::filter(Hugo_Symbol %in% violin_gene)
      violin_mutation_data <- merge(names_mat,violin_mutation_data, by = "Tumor_Sample_Barcode")
      write.csv(violin_mutation_data,paste0(violin_mut_path,"/",violin_gene,"_mutations_in_",disease_filename[j]))
      
      
      
      if (merge_pathways == TRUE) {
        rgenes <- paste0(genes_HGNC,"\\s*?\\.{2}ENSG000",collapse="|")
      }
      else{
        rgenes <- paste0(regulons_violin_gene,"\\s*?\\.{2}ENSG000",collapse="|") 
      }
      
      rgenes_e <- paste0("^(", paste0("CELLLINE|",rgenes), ")")
      
      # #  regression data-set :
      
      mapk_data <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines
                                                    & !(Variant_Classification %in% "Silent"))
      mapk_data <- mapk_data %>% dplyr::filter(Hugo_Symbol %in% violin_gene)
      mapk_data_basis <- mapk_data
      # get whether a mutation is deleterious for the gene function or not so as to approp. change perturbation for CARNIVAL to -1 or 1
      # isDeleterious <-  mapk_data %>% dplyr::select("isDeleterious")
      
      mapk_data <- mapk_data %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol", 
                                               "Variant_Classification","Codon_Change","Protein_Change","isDeleterious")
      
      mapk_data <- add_column(mapk_data, Unique_Mut_ID = 0)
      mapk_data <- mutate(mapk_data, Unique_Mut_ID = sapply(mapk_data$Variant_Classification,
                                                            function(x) grep(x,unique(mapk_data$Variant_Classification))))
      
      colnames(mapk_data) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Codon_Change","Protein_Change","isDeleterious", "Mutation_ID")
      # we dont need detailed mutations in pca matrix as we need to factor for mutation type
      # 
      # 
      names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
      colnames(names_mat) <- c("CELLLINE","CCLE_ID")
      mapk_data <- merge(mapk_data,names_mat, by = "CELLLINE")
      
      
      #####################################################################################################################################################
      # insert start_end positions now in a copy :
      mapk_data_SE <- mapk_data_basis %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol", 
                                                        "Variant_Classification","Codon_Change","Protein_Change","isDeleterious","Start_position","End_position")
      mapk_data_SE <- add_column(mapk_data_SE, Unique_Mut_ID = 0)
      mapk_data_SE <- mutate(mapk_data_SE, Unique_Mut_ID = sapply(mapk_data_SE$Variant_Classification,
                                                                  function(x) grep(x,unique(mapk_data_SE$Variant_Classification))))
      colnames(mapk_data_SE) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Codon_Change","Protein_Change","isDeleterious" ,"Start_position","End_position","Mutation_ID")
      names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
      colnames(names_mat) <- c("CELLLINE","CCLE_ID")
      mapk_data_SE <- merge(mapk_data_SE,names_mat, by = "CELLLINE")
      non_unite_data <- mapk_data_SE 
      #####################################################################################################################################################
      
      cloned <- merge(non_unite_data,filtered_expression_matrix, by = "CELLLINE")
      
      #non_unite_data <- non_unite_data[order(non_unite_data[,'Variant_Classification']), ]
      
      mapk_data <-mapk_data %>% unite(Variant_Classification,Mutation_ID,Protein_Change,CCLE_ID,Variant_Classification, Codon_Change, isDeleterious, sep = "_", remove = TRUE)
      
      colnames(mapk_data) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification")
      
      
      # use expression log2 values for any other multinomial regression to output the data
      r_expr_matrix <- expr_matrix_csv %>% dplyr::select(matches(rgenes_e))
      filtered_expression_matrix_pca <- filtered_expression_matrix_disease  %>% dplyr::select(matches(rgenes_e))
      
      
      # prepare GLM data and run it
      prepare_GLM_data(RNAseq_matrix,
                       rgenes_e,
                       mapk_data,
                       r_expr_matrix,
                       load_other_GLM,
                       GLM_all,
                       reg_type,
                       eXML_dir,
                       key,
                       disease_filename[j],
                       do_GLM,
                       inputs_dir,
                       GLM_predict_user)
      
      
      if (carnival_flag == TRUE || load_other_GLM == TRUE) {
        # prepare Carnival and run it
        if (FEM_user == TRUE) {
          USER_EM <- filtered_expression_matrix_disease
        }
        else{
          USER_EM <- filtered_expression_matrix
        }
        
        prepare_Carnival(mapk_data,
                         USER_EM,
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
                         key,
                         genes_HGNC_bkp,
                         regulons_violin_gene,
                         FC_user,
                         surv,
                         FEM_user,ccle_iterator,condition)
        
      }
      
      mapk_data <- merge(mapk_data, filtered_expression_matrix, by = "CELLLINE", all = TRUE)
      mapk_data$Variant_Classification[is.na(mapk_data$Variant_Classification)] <- "WT"
      mapk_data<-mapk_data[order(mapk_data[,'Variant_Classification']), ]
      setwd(results_dir)
      colnames(mapk_data) <- as.character(colnames(mapk_data))
      colnames(mapk_data) <- sapply(strsplit(colnames(mapk_data), split='..', fixed=TRUE),function(x) (x[1]))
      mapk_data <- mapk_data[,-(1:2),drop=F]
      mapk_data <- t(mapk_data)
      colnames(mapk_data) <- mapk_data[1,]
      mapk_data <- cbind(Gene = rownames(mapk_data), mapk_data)
      rownames(mapk_data) <- NULL
      mapk_data <- mapk_data[-1,,drop=F]
      mapk_data <- mapk_data[,colSums(is.na(mapk_data))<nrow(mapk_data)]
      mapk_data_bkp <- mapk_data
      
      
      # plot the dendrogram:
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
      
      write.csv(t_heat_clust,paste0(results_dir,"/Clustering_data.csv"))
      

      t_heat_clust <- t_heat_clust[ order(row.names(t_heat_clust)), ]  
      #png(file = "CCLE_HEATMAP.png",width = 8, height = 11)
      
      order_rows <- substr(rownames(t_heat_clust),1,1)
      print(order_rows)

      
      print(rownames(t_heat_clust))

      scheme <- c("Missense","Nonsense","In_Frame_Del","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins", "WT")
      
      dend <-  heatmap.2(data.matrix(t_heat_clust), scale = "column",col=bluered(100), breaks=seq(-3, 3, length.out=101),
                         keysize = 1, Rowv=FALSE,
                         trace = "none", density.info = "none", 
                         key.title = "Log2 expression",RowSideColors=as.character(as.numeric(order_rows)),
                         main = "Expression (Log2 Z-transformed) for the regulons VS mutation type",
                         xlab = paste0("Regulons of ", violin_gene, " in ", "CCLE"), 
                         font.lab = 40,ylab = NULL, margins = c(8,20))
      
      legend("topright",      
             legend =  scheme ,
             col = unique(as.numeric(order_rows)), 
             lty= 1,             
             lwd = 10,           
             cex=2.5
      )
      
      print(dend)
      #dev.off()
      
    }
    #******************************************************************************************************************************************
    
    
    setwd(results_dir)
    filtered_mutation_matrix_DL <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines &
                                                                      !(Variant_Classification %in% "Silent") &
                                                                      !(isDeleterious %in% "False"))
    
    filtered_mutation_matrix_DL <- filtered_mutation_matrix_DL %>% dplyr::filter(stringr::str_detect(Hugo_Symbol, exact_genes))
    
    
    grouped_mutations <-  filtered_mutation_matrix %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
    grouped_mutations_DL <-  filtered_mutation_matrix_DL %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
    
    
    
    GM <- grouped_mutations %>% dplyr::summarize(count = n())
    GM_DL <- grouped_mutations_DL %>% dplyr::summarize(count = n())
    
    
    # extract now the pairs of gene-cell line, that is for each gene only the cell line in which it was found mutated
    pairs_GC_CH <- filtered_mutation_matrix %>% dplyr::select(c("Hugo_Symbol","Chromosome", "Tumor_Sample_Barcode",
                                                                "Variant_Classification","Codon_Change","Protein_Change","isDeleterious"))
    
    pairs_GC <- filtered_mutation_matrix %>% dplyr::select(c("Hugo_Symbol","Tumor_Sample_Barcode", "Variant_Classification"))
    ############################################################################################################
    
    filtered_cn_matrix <- cn_csv %>% dplyr::filter(cell_lines %in% filtered_expression_matrix$CELLLINE)
    
    filtered_cn_matrix <- filtered_cn_matrix  %>% dplyr::select(matches(genes_cn))
    mutated_genes <- filtered_mutation_matrix  %>% dplyr::select(matches("Hugo_Symbol"))
    mutated_genes <- unique(mutated_genes[,1])
    
    genes_final <- paste0(mutated_genes,"..ENSG000", collapse="|")
    
    mutated_genes_s <- paste0("^(", paste0("CELLLINE|",genes_final), ")")
    
    #############################################################################################################################
    # calcuate the log fold change of the  median of mutations for each gene across the median on all cell lines in disease)
    
    all_means_mutants <- c()
    all_means_wildtype <- c()
    all_means <- c()
    all_cn_means <- c()
    
    temp <- data.matrix(filtered_expression_matrix[,2:ncol(filtered_expression_matrix)])
    means_E <- apply(temp,2,mean) # applies function 'mean' to 2nd dimension (columns)
    medians_E <- apply(temp,2,median) # applies function 'mean' to 2nd dimension (columns)
    means_medians <- means_E - medians_E
    
    log_fold_change <- c()
    e_val_st_vector <- c()
    
    
    title_cancer = paste0("% done - Analyzing expression profiles in ", disease_filename[j])
    temp_genes <- c() # all the genes appearing mutated but also in the expression matrix
    number_of_mutations <- c() # the number of times a gene was found mutant (maybe in same cell line)
    number_of_unique_mutations <- c() # number of cell lines in which each gene was found mutated (unique)
    
    ################################################ M A I N   L O O P ##################################################
    for (i in 1:length(mutated_genes))        {
      
      lfc <- 0
      means_mutants <- 0
      means_wildtype <- 0
      mutant_expression_cell_lines <- c()
      wildtype_expression_cell_lines <- c()
      temp_column <- c()
      temp_gene0 <-mutated_genes[i]
      
      # make sure a gene found mutated also exists in expression matrix:
      temp_genes_final <- paste0(temp_gene0,"..ENSG000")
      
      temp_mutated_genes_s <- paste0("^",temp_genes_final)
      
      temp_gene0_test <-  filtered_expression_matrix %>% dplyr::select(matches(temp_mutated_genes_s))
      
      if (nrow(temp_gene0_test)  != 0 ) {
        
        temp_gene <- paste0("^(", paste0(temp_gene0, "..ENSG000"), ")")
        
        temp_gene_all <- paste0("^(", paste0("CELLLINE|",temp_gene), ")")
        
        temp_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene_all))
        
        if (temp_gene0 == violin_gene) {
          violin_column <- temp_column[order(temp_column[,'CELLLINE']), ]
        }
        
        if (ncol(temp_column) == 2)
        {
          
          temp_genes <- c(temp_genes, toString(temp_gene0)) # all the genes found mutated but also exist in expression matrix
          mutant_cell_lines0 <-   pairs_GC  %>% dplyr::filter(Hugo_Symbol %in% temp_gene0)
          
          mutant_cell_lines1 <- mutant_cell_lines0 %>% dplyr::select(matches("Tumor_Sample_Barcode", "Variant_Classification"))
          mutant_cell_lines <- as.list(mutant_cell_lines1$Tumor_Sample_Barcode)
          
          number_of_mutations <- c(number_of_mutations,length(mutant_cell_lines0$Tumor_Sample_Barcode))
          number_of_unique_mutations <- c(number_of_unique_mutations,length(unique(mutant_cell_lines0$Tumor_Sample_Barcode)))
          
          mutant_expression_cell_lines <- temp_column  %>% dplyr::filter(CELLLINE %in% mutant_cell_lines0[,2])
          wildtype_expression_cell_lines <- temp_column  %>% dplyr::filter(!(CELLLINE %in% mutant_cell_lines0[,2]))
          
          print("----------------------------------------------------")
          
          # convert to numerical to get means and subtract for log fold change
          mutant_expression_cell_lines <-  as.numeric(na.omit(mutant_expression_cell_lines[,2]))
          wildtype_expression_cell_lines <- as.numeric(na.omit(wildtype_expression_cell_lines[,2]))
          
          means_mutants <- mean(mutant_expression_cell_lines)
          means_wildtype <- mean(wildtype_expression_cell_lines)
          
          
          all_means_mutants = c(all_means_mutants, means_mutants)
          all_means_wildtype = c(all_means_wildtype, means_wildtype)
          
          # this is the expression levels of each gene across all disease-specific  cell lines
          
          gene_expression_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene))
          
          cn_temp_gene <- paste0("^", temp_gene0)
          cn_temp_gene <- paste0(cn_temp_gene,"..{2}")
          gene_cn_column <- filtered_cn_matrix  %>% dplyr::select(matches(cn_temp_gene))
          
          # decide which matrix (Expression or CN) has the less number of disease cell lines to base the test on that subset
          # so as the test can be performed (as equal entries required from both vectors)
          
          l1 <- as.numeric(unlist(gene_expression_column))
          l2 <- as.numeric(unlist(gene_cn_column))
          
          if (ncol(gene_cn_column) == 1) {
            
            all_means <- c(all_means,median(l1))
            all_cn_means <- c(all_cn_means,median(l2))
          }
          
          ##############################################################################################################
          
          #}
        }
      }
    } # end of mutated genes iterator loop
    
    ttest_matrix = cbind(all_means_mutants,all_means_wildtype)
    p_val_st <- 0 
    p_val <- 0 
    e_val_st <- 0
    pval <- 0
    mean_of_differences <- 0
    
    tt <- NULL
    
    tt <- t.test(ttest_matrix[,1], ttest_matrix[,2], paired = TRUE)
    p_val <- round(tt$p.value, digits = 5)
    mean_of_differences <- round(tt$estimate, digits = 5)
    
    
    #print("Spearman test between the median of a Gene expression against the mean of its CN across all disease cell lines:")
    stest_matrix = cbind(all_means,all_cn_means)
    
    ST <- cor.test(stest_matrix[,1],stest_matrix[,2] , method = "spearman", exact=FALSE )
    
    p_val_st <- round(ST$p.value, digits = 5)
    e_val_st <- round(ST$estimate, digits = 5)
    
    coverage = round((length(filtered_expression_matrix$CELLLINE)/length(disease_cell_lines))*100, digits = 1)
    
    amplified_genes = length(all_cn_means)
    
    log2_fold_change = ttest_matrix[,1] -  ttest_matrix[,2]
    
    
    main1 = paste0("(Log2) Fold change (MT-WL) of MAPK/ERK pathway genes in found mutant VS WT cell lines in ",disease_filename[j],
                   ".\n ----------------------------------Statistics---------------------------------------- \n | Genes in Pathway: ", length(genes_HGNC), paste0(". Regulons of ", violin_gene, ":"), 
                   length(regulons_violin_gene),"."
                   ,     "\n | P-value of pairwise T-test: ", p_val, ". Mean of differences : ",
                   mean_of_differences , "\n | P-value of Spearman test (expression VS CN): ",p_val_st ,
                   ". Estimate of Spearman correlation (rho): ",  e_val_st, "\n | Number of Cell-Lines found: ",length(filtered_expression_matrix$CELLLINE) ,
                   "/", length(disease_cell_lines) , ". Lineage Coverage : ",  coverage, "%", " .Genes found mutated/amplified: ", length(all_means_mutants) ,"/",
                   amplified_genes, ".",
                   "\n | (source data-sets: DepMap Public ", dataset_version, ")", "\n ---------------------------------------------------------------------------------")
    
    
    gene_labels <- paste0(as.list(temp_genes), "(", as.list(number_of_unique_mutations),"/",length(filtered_expression_matrix$CELLLINE),")")
    
    log2_fold_change_df <- data.frame(temp_genes,number_of_mutations,number_of_unique_mutations,
                                      ttest_matrix[,1],ttest_matrix[,2],log2_fold_change,gene_labels)
    names(log2_fold_change_df) <- c("Gene", "Number_of_Mutations", "Number_of_Unique_Mutations", "Mutant",
                                    "WT", "Log2_Fold_difference", "LABELS")
    
    log2_fold_change_df <- log2_fold_change_df %>% arrange(desc(number_of_mutations))
    
    Number_of_Mutations <- log2_fold_change_df$Number_of_Mutations
    
    ########################################### MUTATIONS OF ALL GENES GRAPH and T-test #################################################
    # just get the labels as : name of gene + (number_of_unique_mutations)
    
    
    sp2 <- ggplot(log2_fold_change_df,  aes(x=factor(Gene, levels=Gene), y=Log2_Fold_difference,label = LABELS,size = Number_of_Mutations)) +
      
      geom_point(color = dplyr::case_when(log2_fold_change_df$Log2_Fold_difference > 1 ~ "#FF0000",
                                          log2_fold_change_df$Log2_Fold_difference < -1 ~ "#FF0000",
                                          TRUE ~ "#00CC00"), alpha = 0.8) +
      geom_hline(yintercept = 1,linetype="dotted") +
      geom_hline(yintercept = -1,linetype="dotted") +
      
      geom_hline(yintercept = 2,linetype="dashed") +
      geom_hline(yintercept = -2,linetype="dashed") +
      geom_text_repel(label = ifelse(log2_fold_change_df$Number_of_Mutations >= 1 &
                                       log2_fold_change_df$Log2_Fold_difference < 1 &
                                       log2_fold_change_df$Log2_Fold_difference > -1, as.character(log2_fold_change_df$LABELS) , "" ), size = 2) +
      geom_text_repel( data          = subset(log2_fold_change_df, Log2_Fold_difference > 1),
                       nudge_y       = 16 - subset(log2_fold_change_df, Log2_Fold_difference > 1)$Log2_Fold_difference,
                       size          = 2,
                       box.padding   = 1.5,
                       point.padding = 0.5,
                       force         = 0.5,
                       segment.size  = 0.5,
                       segment.color = "grey50",
                       direction     = "y") +
      geom_label_repel(data         = subset(log2_fold_change_df, Log2_Fold_difference < -1),
                       nudge_y       = -16 - subset(log2_fold_change_df, Log2_Fold_difference < -1)$Log2_Fold_difference,
                       size          = 2,
                       box.padding   = 0.5,
                       point.padding = 0.5,
                       force         = 0.5,
                       segment.size  = 0.5,
                       segment.color = "grey50",
                       direction     = "y") +
      
      scale_x_discrete(expand = expand_scale(mult = c(0.005, .05))) +
      scale_y_continuous(expand = expand_scale(mult = c(0.005, .01)))  +
      theme(axis.text.x=element_text(size=5, angle=90,hjust=0.95,vjust=0.2),plot.title = element_text(size = 9)) + ggtitle(main1)
    
    
    
    ######################################### MUTATION TYPES GRAPH ###############################################
    
    main2 = paste0(" Mutation variation per gene in ",disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    
    GM_df <- data.frame(GM$Hugo_Symbol,GM$Variant_Classification,GM$count)
    
    names(GM_df) <- c("Hugo_Symbol", "Variant_Classification", "count")
    
    
    
    sp3<-  ggplot(GM_df,aes(x = Hugo_Symbol, y = count, fill = Variant_Classification)) + geom_bar(stat = "identity") + ggtitle(main2) +
      
      theme(axis.text.x=element_text(size=6, angle=90,hjust=0.95,vjust=0.2)) +
      geom_text(aes(label=count),size = 3, hjust = 0.5, vjust = 3, position ="stack")
    
    #################################### FILTERED ON IS_DELETERIOUS OR NOT #####################################
    
    main2 = paste0(" Deleterious Mutation variation per gene in ",disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    
    GM_df_DL <- data.frame(GM_DL$Hugo_Symbol,GM_DL$Variant_Classification,GM_DL$count)
    
    names(GM_df_DL) <- c("Hugo_Symbol", "Variant_Classification", "count")
    
    sp4<-  ggplot(GM_df_DL,aes(x = Hugo_Symbol, y = count, fill = Variant_Classification)) + geom_bar(stat = "identity") + ggtitle(main2) +
      
      theme(axis.text.x=element_text(size=6, angle=90,hjust=0.95,vjust=0.2)) +
      geom_text(aes(label=count),size = 3, hjust = 0.5, vjust = 3, position ="stack")
    
    
    #################################### Mean of Expression - Median of Expression #####################################
    
    main2 = paste0("Difference of mean and median values of expression per gene in ",
                   disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    GM_MM <- data.frame(names(filtered_expression_matrix[,2:ncol(filtered_expression_matrix)]),means_medians)
    
    
    
    colnames(GM_MM) <- c("Gene_ID", "Difference")
    
    GM_MM$Gene_ID <- as.numeric(factor(GM_MM$Gene_ID))
    
    sp5<-  ggplot(GM_MM,aes(x = Gene_ID, y = Difference) ) +
      geom_point(color="blue", size=2)+
      geom_smooth(method = "lm", se = TRUE) +
      ggtitle(main2) +
      theme(axis.text.x=element_text(size=6, angle=90,hjust=0.95,vjust=0.3))
    
    pairs_GC_alp <-   pairs_GC_CH[order(pairs_GC[,'Tumor_Sample_Barcode']), ]
    
    # for violin_gene
    pairs_GC_alp_bkp <-  pairs_GC_alp
    pairs_GC_alp <- pairs_GC_alp %>% dplyr::filter(Hugo_Symbol %in% violin_gene & Tumor_Sample_Barcode %in% violin_column$CELLLINE)
    
    
    pairs_GC_alp <- add_column(pairs_GC_alp, Expression = 0)
    pairs_GC_CH <- add_column(pairs_GC_CH, Expression = 0)
    
    pairs_GC_alp  <- pairs_GC_alp   %>% dplyr::mutate(Hugo_Symbol = disease_filename[j])
    
    # now pairs_GC_alph has all cellines, mutation types and corresponding expression levels for the violin gene
    
    if (!is.null(violin_column)) {
      
      violin_column <- violin_column %>% dplyr::filter(CELLLINE %in% pairs_GC_alp$Tumor_Sample_Barcode)
      # for-loop to fill-in the expression levels of all mutated cell lines for the "Violin-Gene" (usually TP53)
      j2 = 1
      for (i in 1: length(pairs_GC_alp$Tumor_Sample_Barcode)) {
        
        index = match(pairs_GC_alp$Tumor_Sample_Barcode[i],violin_column$CELLLINE)
        pairs_GC_alp$Expression[j2] = violin_column[index,2]
        
        j2 = j2 + 1
      }
    }
    # same but for all genes
    j2 = 1
    
    print("Analyzing mutation profiles...")
    
    for (i in 1: length(pairs_GC_CH$Tumor_Sample_Barcode)) {
      temp_gene <- pairs_GC_CH$Hugo_Symbol[i]
      
      temp_gene_all <- paste0("^(", paste0("CELLLINE|",temp_gene), ")")
      
      
      temp_column <- NULL
      temp_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene_all))
      if (ncol(temp_column) == 2)     {
        
        index = match(pairs_GC_CH$Tumor_Sample_Barcode[i],temp_column$CELLLINE)
        
        pairs_GC_CH$Expression[j2] = temp_column[index,2]
      }
      else{
        pairs_GC_CH$Expression[j2] = NA
        
      }
      
      j2 = j2 + 1
    }

    ###### (Box - Violin - Scatter Plot) on all mutated genes across cell lines and types of mutations versus expression level ##########
    dodge <- position_dodge(width = 0.4)
    
    
    dodge <- position_dodge(width = 0.4)
    
    main6_1 = paste0("Violin and boxplots of types of mutations across all mutated genes found in ",
                     disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    main6_2 = paste0("Violin and boxplots of types of mutations (annotated with Gene Hugo Symbols) across all mutated genes found in ",
                     disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    main6_3 = paste0("Violin and boxplots of types of mutations (annotated with cell line names) across all mutated genes found in ",
                     disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
    
    sp6_1 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
      #scale_fill_viridis_d( option = "D")+
      
      geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
      geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +
      
      
      #geom_label_repel(aes(label=total$cell_lines),
      #                 box.padding   = 0.5,
      #                 point.padding = 0.005, size = 1.8) +
      
      
      geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
      ylab(  c("Expression (log2 values)")  )  +
      xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +
      
      font("xylab",size=20)+
      font("xy",size=20)+
      font("xy.text", size = 20) +
      font("legend.text",size = 20) +
      theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02)) +
      ggtitle(main6_1)
    
    sp6_2 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
      #scale_fill_viridis_d( option = "D")+
      
      geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
      geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +
      geom_label_repel(aes(label=pairs_GC_CH$Hugo_Symbol),
                       box.padding   = 0.5,
                       point.padding = 0.005,
                       segment.color = 'grey50', size = 1.5) +
      geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
      ylab(  c("Expression (log2 values)")  )  +
      xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +
      
      font("xylab",size=20)+
      font("xy",size=20)+
      font("xy.text", size = 20) +
      font("legend.text",size = 20) +
      theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02)) +
      ggtitle(main6_2)
    
    
    
    sp6_3 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
      #scale_fill_viridis_d( option = "D")+
      
      geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
      geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +
      geom_label_repel(aes(label=pairs_GC_CH$Tumor_Sample),
                       box.padding   = 0.5,
                       point.padding = 0.005,
                       segment.color = 'grey50', size = 1) +
      geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
      ylab(  c("Expression (log2 values)")  )  +
      xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +
      
      font("xylab",size=20)+
      font("xy",size=20)+
      font("xy.text", size = 20) +
      font("legend.text",size = 20) +
      theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
      ggtitle(main6_3)
    
    
    
    ###### (Box - Violin - Scatter Plot) on violin_gene across cell lines and types of mutations versus expression level ##########
    genes <- paste0(violin_gene,"\\s*?\\.{2}ENSG000",collapse="|")
    
    genes0 <- paste0("\\b",violin_gene,"\\s*?\\.{2}\\b", collapse="|")
    genes_e <- paste0("^(", paste0("CELLLINE|",genes), ")")
    genes_cn <- paste0("^(", paste0("cell_lines|",genes0), ")")
    
    test_data <- data.frame()
    
    
    # expression
    test_data <- expr_matrix_csv  %>% dplyr::select(matches(genes_e))
    test_data$WT <- ifelse(expr_matrix_csv$CELLLINE %in% pairs_GC_alp$Tumor_Sample_Barcode, "MT", "WT")
    colnames( test_data) <- c("cell_lines", "Expression_log2", "State")
  
    #print(test_data)
    # copy number
    cn_gene <- data.frame()
    cn_gene <- cn_csv  %>% dplyr::select(matches(genes_cn))
    
    if (ncol(cn_gene) == 1) # no CN data was found for the violin gene
    {
      cn_gene[,2] <- rep(0,nrow(cn_gene)) # filling then CN data with 0 zeroes
    }
    else
    {
      
      # plot CN across Expression levels for violin gene only
      violin_gene_CN <- NULL
      violin_gene_E <- NULL
      
      violin_gene_CN <- cn_gene
      
      violin_gene_E <- expr_matrix_csv %>% dplyr::select(matches(genes_e))
      
      if (ncol(violin_gene_E)==2 & ncol(violin_gene_CN)==2) {
        colnames(violin_gene_E) <- c("CELLLINE",violin_gene)
        colnames(violin_gene_CN) <- c("CELLLINE",violin_gene)
        
        #violin_gene_CN_E <- NULL
        
        violin_gene_CN_E <- merge(violin_gene_CN,violin_gene_E, by = "CELLLINE")
       # print(violin_gene_CN_E)
        colnames(violin_gene_CN_E) <- c("CELLLINE","CN_log2","Expression_log2")
        violin_gene_CN_E <- violin_gene_CN_E %>% dplyr::select("CN_log2","Expression_log2")
        
        
        violin_gene_CN_E_plot <- ggplot(violin_gene_CN_E, aes(x=CN_log2, y=Expression_log2))+
          scale_x_continuous(breaks=seq(0,80,1)) +
          geom_point()+
          ggtitle(paste0("Copy number (log2) versus expression (log2) for ", violin_gene," in all cell lines"))
        print(violin_gene_CN_E_plot)
      }
      
    }
    
    cn_gene <-  cn_gene  %>% dplyr::select(c(1,2))
    colnames(cn_gene) <- c("cell_lines", "CN")
    cn_gene$CN_S <- ifelse(cn_gene$CN >= 1, "Amplification", "Deletion")
    
    # merging dataframes to one
    total_o <- merge(test_data, cn_gene, by = "cell_lines")
    
    pairs_GC_alp_temp <- NULL
    pairs_GC_alp_temp <- pairs_GC_alp
    colnames(pairs_GC_alp_temp)[3] <- "cell_lines"
    
    pairs_GC_alp_temp<- pairs_GC_alp_temp %>% dplyr::select("cell_lines","Variant_Classification","isDeleterious", "Protein_Change")
    
    
   
    total <- merge(total_o, pairs_GC_alp_temp, by = "cell_lines", all = TRUE)
    
    total <- total[!is.na(total$State),]
    
    total  <- mutate(total, Variant_Classification = ifelse(is.na(Variant_Classification), "WT", as.character(Variant_Classification)))
    write.csv(calls,"calls.csv")

    total <- merge(total,calls, by = "cell_lines")
    total  <- mutate(total, isHotspot = ifelse(Protein_Change %in% hotspots_bkp, TRUE, FALSE))
    write.csv(total,"total.csv")
    
    
    
    if (!is.null(violin_column)) {
      ####### is Hotspot ######
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                     " \n (source data-sets: DepMap Public  ", dataset_version, ")")
      
      total2 <- total[!is.na(total$GYSTIC_CALLS),] 
      
      sp777 <- ggplot(data = total2,aes(x = State, y = Expression_log2, fill = isHotspot))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) +
        stat_compare_means(method = "anova", label.y = 20, size = 5)  + 
        stat_compare_means(label.y = 21, size = 5) +
        stat_compare_means(method = "t.test",label.y = 22, size = 5) +
        stat_compare_means(method = "wilcox.test",label.y = 23, size = 5)
        
      print(sp777)
      plot_list <- c(plot_list,sp777)
      ggsave(filename="CCLE_WT_MT_ishotspot.png", plot=sp777)
      ####### is deleterious ######
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                     " \n (source data-sets: DepMap Public  ", dataset_version, ")")
      

      
      sp777 <- ggplot(data = total2,aes(x = State, y = Expression_log2, fill = isDeleterious))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) +
        stat_compare_means(method = "anova", label.y = 20, size = 5)  + 
        stat_compare_means(label.y = 21, size = 5) +
        stat_compare_means(method = "t.test",label.y = 22, size = 5) +
        stat_compare_means(method = "wilcox.test",label.y = 23, size = 5)
      
      print(sp777)
      plot_list <- c(plot_list,sp777)
      ggsave(filename="CCLE_WT_MT_isdel.png", plot=sp777)
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                     " \n (source data-sets: DepMap Public  ", dataset_version, ")")
      
      total2 <- total[!is.na(total$GYSTIC_CALLS),] 
      #ggplot(aes(x=reorder(carrier,-speed, na.rm = TRUE)
      # ggplot(data = total2,aes(x = State, y = Expression_log2,      
 
      sp777 <- ggplot(data = total2,aes(x = State,y = Expression_log2, fill = GYSTIC_CALLS))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) +
        stat_compare_means(method = "anova", label.y = 12, size = 8)  + 
        stat_compare_means(label.y = 15, size = 8)+
        #stat_compare_means(method = "t.test",label.y = 22, size = 8) +
        #stat_compare_means(method = "wilcox.test",label.y = 16, size = 8) +
        stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")
      
      
      print(sp777)
      plot_list <- c(plot_list,sp777)
      ggsave(filename="CCLE_WT_MT_gystic.png", plot=sp777)
      ####### is deleterious ######
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                     " \n (source data-sets: DepMap Public  ", dataset_version, ")")
      
      
      
      
      
    }
    if (!is.null(violin_column)) {
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels and types for ", violin_gene, " in all cell lines",
                     " \n (source data-sets: DepMap Public  ", dataset_version, ")")
      
      
      sp77 <- ggplot(data = total,aes(x = Variant_Classification, y = Expression_log2, fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        #geom_point(shape = ifelse(total$CN_S == "Amplification", 8, 21),size=ifelse(total$CN_S == "Amplification", 3, 1),
        #           position = dodge, color= ifelse(total$CN_S == "Amplification", "red", "yellow"),alpha=1) +
        #geom_label_repel(aes(label=total$cell_lines),
        #                 box.padding   = 0.5,
        #                 point.padding = 0.005, size = 1.8) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7) +
        stat_compare_means(method = "anova", label.y = 20, size = 5)  + 
        stat_compare_means(label.y = 21, size = 5)
      
      print(sp77)
      plot_list <- c(plot_list,sp77)
      ggsave(filename="CCLE.png", plot=sp77)
    }
    
    if (!is.null(violin_column)) {
      
      dodge <- position_dodge(width = .4)
      main7 = paste0("Violin and boxplots of types of mutations for ", violin_gene, " in ",
                     disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
      
      
      colnames(names_mat)[1] <- "Tumor_Sample_Barcode"
      pairs_GC_alp <- merge(names_mat, pairs_GC_alp, by = "Tumor_Sample_Barcode")
      
      
      sp7 <- ggplot(data = pairs_GC_alp,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+
        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
        geom_label_repel(aes(label=pairs_GC_alp$CCLE_ID),
                         box.padding   = 0.5,
                         point.padding = 0.005, size = 3) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0(violin_gene," Variant Classification in ", disease_filename[j])) ) +
        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main7)
      print(sp7)
      plot_list <- c(plot_list,sp7)
      
    }
    all_cancers <- rbind(all_cancers,pairs_GC_alp)
    
    
    #pairs_GC_CH_bkp <- NULL
    pairs_GC_CH_bkp <- pairs_GC_CH
    
    # map broad ids to cell line ids
    colnames(names_mat)[1] <- "Tumor_Sample_Barcode"
    pairs_GC_CH_bkp <- merge(names_mat, pairs_GC_CH_bkp, by = "Tumor_Sample_Barcode")
    pairs_GC_CH_bkp <- pairs_GC_CH_bkp %>% dplyr::select(c("Hugo_Symbol","Chromosome", "CCLE_ID", 
                                                           "Variant_Classification", "Codon_Change" , "Protein_Change" , "Expression"))
    colnames(pairs_GC_CH_bkp)[3] <- "Tumor_Sample_Barcode"
    
    all_cancers_all_genes <- rbind(all_cancers_all_genes,pairs_GC_CH)
    
    print(sp2)
    print(sp3)
    print(sp4)
    print(sp5)
    print(sp6_1)
    plot_list <- c(plot_list,sp2,sp3,sp4,sp5,sp6_1)
    
    pairs_GC_CH  <-  pairs_GC_CH %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
    pairs_GC_CH <- pairs_GC_CH %>% dplyr::summarize(count = n())
    
    pairs_GC_CH <- pairs_GC_CH  %>% spread(key=Variant_Classification, value=count)
    #pairs_GC_CH[is.na(pairs_GC_CH)] <- 0
    pairs_GC_CH$Hugo_Symbol <- NULL
    pairs_GC_CH <- pairs_GC_CH[2:nrow(pairs_GC_CH),]
    
    # clear zero columns:
    #
    # if (nrow(pairs_GC_CH) <= 1) {
    #   print("No mutations found for the signature genes. Cannot do mutation corrplot.")
    # }
    # else{
    #   pairs_GC_CH <- pairs_GC_CH[,which(colSums(pairs_GC_CH)!=0)]
    #   title2 <- paste0("Correlation matrix of mutation variation for all genes across ",
    #                    disease_filename[j], " (source data-sets: DepMap Public  ", dataset_version, ")")
    #   
    #   
    #   ACCM2 <-cor( data.matrix(pairs_GC_CH), method = "spearman")
    #   write.csv(pairs_GC_CH,paste0(results_dir,"/pairs_GC_CH.csv"))
    #   sp10 <- corrplot(ACCM2,type = "lower", col = cm.colors(100), title = title2, mar=c(0,0,1,0))
    #   print(sp10)
    #   plot_list <- c(plot_list,sp10)
    #   
    # }
    # if (disease_name != "all_cell_lines") {
    #   # dynamic_graphs()
    #   dynamic_graphs_res <-dynamic_graphs(pairs_GC_CH_bkp,
    #                                       violin_column,
    #                                       violin_gene,
    #                                       regulons_violin_gene,
    #                                       dynamic_graphs_dir,
    #                                       disease_filename[j],
    #                                       skip_C)
    #   
    #   print(dynamic_graphs_res$sp_g)
    # }
    ###################################### IGRAPH PLOT ####################################
    
    # now plot cell lines heatmap versus mutations of TP53, the expression of which is the color
    if (!is.null(violin_column)) {
      
      sp12 <- ggplot(data = pairs_GC_alp, aes(x = Variant_Classification, y = CCLE_ID)) +
        geom_tile(aes(fill = Expression)) +
        geom_text(aes(label=paste0(Protein_Change,"\n",Codon_Change)),size = 2) +
        theme(axis.text.x=element_text(size=15, angle=90,hjust=0.95,vjust=0.02))+
        scale_fill_distiller(palette = "Spectral") +
        ggtitle(paste0("Heatmap of expression levels of ", violin_gene,
                       " across all cell lines and mutation variants in ", disease_filename[j],
                       "\n (source data-sets: DepMap Public  ", dataset_version, ")"))
      
      print(sp12)
      plot_list <- c(plot_list,sp12)
    }
    
    sp13 <- ggplot(data = pairs_GC_CH_bkp, aes(x = Variant_Classification, y = Tumor_Sample_Barcode)) +
      geom_tile(aes(fill = Expression)) +
      theme(axis.text.x=element_text(size=15, angle=90,hjust=0.95,vjust=0.02))+
      scale_fill_distiller(palette = "Spectral") +
      ggtitle(paste0("Heatmap of expression levels of all mutated genes across all cell lines and mutation variants in ",
                     disease_filename[j], "\n (source data-sets: DepMap Public  ", dataset_version, ")"))
    
    print(sp13)
    plot_list <- c(plot_list,sp13)
    
    
    dev.off()
  } # end of cancer type iterator loop
  
  
  if (length(disease_filename) > 1000) {
    pdffile = paste0("pancancer_",dataset_version,".pdf")
    pdf(file = pdffile, height = 10, width = 18, onefile = TRUE)
    
    pancancer_graphs_res <- pancancer_graphs(violin_gene,
                                             all_cancers_all_genes,
                                             all_cancers,
                                             pairs_GC_CH_bkp,
                                             skip_C,
                                             results_dir,
                                             wd,
                                             mut_matrix_csv,
                                             expr_matrix_csv,
                                             mapk_data_basis,
                                             filtered_expression_matrix,
                                             pca_objects,
                                             dataset_version)
    
    print(pancancer_graphs_res$pca_whole_res$pca_plot_default)
    print(pancancer_graphs_res$sp9)
    print(pancancer_graphs_res$sp11)
    print(pancancer_graphs_res$sp12_1)
    print(pancancer_graphs_res$sp_violin)
    
    dev.off()
  }

  
  print("Finished scan.")
  end <- Sys.time()
  time <- end - start
  print(time)
  print("---------------END OF ANALYSIS--------------")
  
  # end of function
}
