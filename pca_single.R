# pca per cancer module

###########################################################################
# create the PCA matrix:

pca_single <- function(violin_gene,
                       disease_filename_j,
                       mapk_data_basis,
                       filtered_expression_matrix,
                       results_dir,
                       pca_objects) {
  
  
  
  filtered_expression_matrix <- filtered_expression_matrix %>% dplyr::select(-matches(violin_gene))
  
  
  print(paste0("Calculating PCA on ", violin_gene," mutation status in ", disease_filename_j))
  mapk_data_pca <- mapk_data_basis %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  colnames(mapk_data_pca) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification")
  
  
  mapk_data_pca <- unique(mapk_data_pca)
  r_expr_matrix_pca <- unique(filtered_expression_matrix)
  pca_matrix <- merge(mapk_data_pca, r_expr_matrix_pca, by = "CELLLINE", all = TRUE)
  
  # Missense_Mutation
  # pca_matrix <- pca_matrix %>%  dplyr::filter(Variant_Classification %in% "Missense_Mutation")
  
  
  
  pca_matrix<-pca_matrix[order(pca_matrix[,'Variant_Classification']), ]
  pca_matrix$Variant_Classification[is.na(pca_matrix$Variant_Classification)] <- "WT"
  pca_matrix <- mutate(pca_matrix, Variant_Classification = as.character(Variant_Classification))
  

  pca_matrix[is.na(pca_matrix)] <- 0
  
  
  colnames(pca_matrix) <- sapply(strsplit(colnames(pca_matrix), split='..', fixed=TRUE),function(x) (x[1]))
  meta_pca <- as.data.frame(pca_matrix)
  meta_pca$Variant_Classification[pca_matrix$Variant_Classification == 0] <- "WT"
  pca_matrix <- meta_pca[,-(2:3),drop=F]
  pca_matrix <- as.data.frame(data.matrix(pca_matrix))
  pca_matrix <- pca_matrix[,which(colSums(pca_matrix) != 0)]
  pca_matrix <- pca_matrix[,-(1),drop=F]
  
  write.csv(meta_pca,paste0(results_dir,"/meta.csv"))
  write.csv(meta_pca,paste0(results_dir,"/pca_matrix.csv"))

  pca_object <- prcomp(pca_matrix,scale=TRUE)

  # now do the plots
  pca_plot <- plot(ggplot2::autoplot(pca_object, data = meta_pca,loadings = TRUE, loadings.colour = 'blue',
                                loadings.label = TRUE, colour = "Variant_Classification"))
  pca_plot_frame <- plot(ggplot2::autoplot(pca_object, data = meta_pca, shape= "Variant_Classification", frame=T, colour = "Variant_Classification"))
  pca_plot_default <- plot(ggplot2::autoplot(pca_object, data = meta_pca, colour = "Variant_Classification"))

  # save pca object
  saveRDS(pca_object, file = paste0(pca_objects,"/",disease_filename_j,"_pca.rds"))

  return_list <- list(pca_plot,pca_plot_frame, pca_plot_default)
  return(return_list)
}
############################################################################