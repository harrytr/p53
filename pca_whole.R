# pca_whole module

pca_whole <- function(violin_gene,
                      mut_matrix_csv,
                      expr_matrix_csv,
                      mapk_data_basis,
                      filtered_expression_matrix,
                      results_dir,
                      pca_objects) {
  
  filtered_expression_matrix <- filtered_expression_matrix %>% dplyr::select(-matches(violin_gene))
  
  
  print(paste0("Calculating PCA of ", violin_gene, " mutation status in whole transcriptome..."))
  # now include all samples in PCA and all genes
  mapk_data_basis <- mut_matrix_csv %>% dplyr::filter(!(Variant_Classification %in% "Silent"))
  mapk_data_basis <- mapk_data_basis %>% dplyr::filter(Hugo_Symbol %in% violin_gene)
  
  
  mapk_data_basis <- unique(mapk_data_basis)
  
  expr_matrix_csv_temp <- unique(expr_matrix_csv)
  
  mapk_data_pca <- mapk_data_basis %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  colnames(mapk_data_pca) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification")
  
  pca_matrix <- merge(mapk_data_pca, expr_matrix_csv_temp, by = "CELLLINE", all = TRUE)
  pca_matrix<-pca_matrix[order(pca_matrix[,'Variant_Classification']), ]
  pca_matrix$Variant_Classification[is.na(pca_matrix$Variant_Classification)] <- "WT"
  pca_matrix <- mutate(pca_matrix, Variant_Classification = as.character(Variant_Classification))
  pca_matrix[is.na(pca_matrix)] <- 0
  
  #colnames(pca_matrix) <- as.character(colnames(pca_matrix))
  colnames(pca_matrix) <- sapply(strsplit(colnames(pca_matrix), split='..', fixed=TRUE),function(x) (x[1]))
  
  meta_pca <- as.data.frame(pca_matrix)
  meta_pca$Variant_Classification[pca_matrix$Variant_Classification == 0] <- "WT"
  pca_matrix <- meta_pca[,-(2:3),drop=F]
  pca_matrix <- as.data.frame(data.matrix(pca_matrix))
  pca_matrix <- pca_matrix[,which(colSums(pca_matrix) != 0)]
  
  pca_matrix <- pca_matrix[,-(1),drop=F]
  write.csv(meta_pca,paste0(results_dir,"/whole_pca_matrix.csv"))
  
  pca_object <- prcomp(pca_matrix,scale=TRUE)
  
  # now do the plots
  #pca_plot <- ggplot2::autoplot(pca_object, data = meta_pca,loadings = TRUE, loadings.colour = 'blue',
  # loadings.label = TRUE, colour = "Variant_Classification")
  #pca_plot_frame <- ggplot2::autoplot(pca_object, data = meta_pca, shape= "Variant_Classification", frame=T, colour = "Variant_Classification")
  
  pca_plot_default <- plot(ggplot2::autoplot(pca_object, data = meta_pca, colour = "Variant_Classification") + theme(text = element_text(size = 20)))
  
  print(pca_plot_default)
  #print(pca_plot_frame)
  saveRDS(pca_object, file = paste0(pca_objects,"/_whole_transcriptome_pca.rds"))
  return_list <- list(pca_plot_default)
  return(return_list)
  
  
}