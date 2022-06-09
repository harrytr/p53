

# pancancer graphs module 

pancancer_graphs <- function(violin_gene,
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
                             dataset_version) {
  
  source(paste0(wd,'//pca_whole.R'))
  
  
  all_cancers_all_genes_bkp <- NULL
  all_cancers_all_genes_bkp <- all_cancers_all_genes
  all_cancers_all_genes_bkp <- all_cancers_all_genes_bkp %>% dplyr::filter(Hugo_Symbol %in% violin_gene)
  
  all_cancers_f  <-  all_cancers %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
  all_cancers_f <- all_cancers_f %>% dplyr::summarize(count = n())
  
  all_cancers_all_genes <- all_cancers_all_genes  %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
  all_cancers_all_genes <- all_cancers_all_genes  %>% dplyr::summarize(count = n())
  
  all_cancers_f2 <- all_cancers_f  %>% spread(key=Variant_Classification, value=count)
  all_cancers_f2[is.na(all_cancers_f2)] <- 0
  all_cancers_f2$Hugo_Symbol <- NULL
  all_cancers_f2 <- all_cancers_f2[2:nrow(all_cancers_f2),]
  
  
  all_cancers_all_genes <-  all_cancers_all_genes %>% spread(key=Variant_Classification, value=count)
  all_cancers_all_genes[is.na(all_cancers_all_genes)] <- 0
  all_cancers_all_genes$Hugo_Symbol <- NULL
  all_cancers_all_genes <- all_cancers_all_genes[2:nrow(all_cancers_all_genes),]
  
  
  pairs_GC_CH_bkp  <-  pairs_GC_CH_bkp %>% dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol)
  pairs_GC_CH_bkp <- pairs_GC_CH_bkp %>% dplyr::summarize(count = n())
  
  pairs_GC_CH_bkp <-pairs_GC_CH_bkp  %>% spread(key=Hugo_Symbol, value=count)
  pairs_GC_CH_bkp[is.na(pairs_GC_CH_bkp)] <- 0
  print(pairs_GC_CH_bkp)
  
  write.csv(pairs_GC_CH_bkp,paste0(results_dir,"/pairs_GC_CH_bkp.csv"))
  
  ACCM <-cor( data.matrix(pairs_GC_CH_bkp), method = "spearman")
  
  
  # if (skip_C == FALSE) {
  #   
  #   pca_whole_res <- pca_whole(violin_gene,
  #                         mut_matrix_csv,
  #                         expr_matrix_csv,
  #                         mapk_data_basis,
  #                         filtered_expression_matrix,
  #                         results_dir,
  #                         pca_objects)
  #     
  # }
  
  ############################################################################
  
  ###### separating title page for each cancer type #####
  a = paste0("Pan-cancer results for ",violin_gene , " and all other genes")
  
  
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(1,4,a, pos=4)
  
  
  #####################################################
  title1 <- paste0("Correlation matrix of TP53 mutation variation across all test cancer types (source data-sets: DepMap Public  ", dataset_version, ")")
  title3 <- paste0("Correlation matrix of mutation variation across all genes and cancer types (source data-sets: DepMap Public  ", dataset_version, ")")
  
  ACCM <-cor( data.matrix(all_cancers_f2), method = "spearman")
  sp9 <- corrplot(ACCM, method="color", type = "lower", title = title1,mar=c(0,0,1,0))
  
  
  ACCM3 <-cor( data.matrix(all_cancers_all_genes), method = "spearman")
  sp11 <- corrplot(ACCM3, method="color", type = "lower", , title = title3,mar=c(0,0,1,0))
  
  
  sp12_1 <- ggplot(data =all_cancers_all_genes_bkp, aes(x = Variant_Classification, y = Tumor_Sample_Barcode)) +
    geom_tile(aes(fill = Expression)) +
    theme(axis.text.x=element_text(size=15, angle=90,hjust=0.95,vjust=0.02))+
    ggtitle(paste0("Heatmap of expression levels of ", violin_gene, " across all cell lines and mutation variants in all cancers",
                   "\n (source data-sets: DepMap Public  ", dataset_version, ")"))
  
  
  ############################# violin and boxplots of violin_Gene across all run cancer types ###########
  
  dodge <- position_dodge(width = .4)
  main7 = paste0("Violin and boxplots of expression levels for ", violin_gene, " in all cancer types",
                 " \n (source data-sets: DepMap Public  ", dataset_version, ")")
  
  
  sp_violin <- ggplot(data = all_cancers,aes(x = tools::file_path_sans_ext(Hugo_Symbol), y = Expression, fill = Variant_Classification))+
    #scale_fill_viridis_d( option = "D")+
    geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
    geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
    geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
    ylab(  c("Expression (log2 values)")  )  +
    xlab(  c(paste0(violin_gene," in all cancer types tested") )) +
    font("xylab",size=10)+
    font("xy",size=10)+
    font("xy.text", size = 10) +
    font("legend.text",size = 10) +
    theme(axis.text.x=element_text(size=15, angle=90,hjust=0.95,vjust=0.02))+
    ggtitle(main7)
  
  
  return_list <- list(pca_whole_res$pca_plot_default,plot(sp9),plot(sp11),plot(sp12_1),plot(sp_violin))
  return(return_list)
  
}