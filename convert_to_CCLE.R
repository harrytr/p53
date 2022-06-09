convert_to_CCLE <- function()
  
{
  print("Loading the data...")
  #load ("TCGA_RNA.RData")
  library(stringr)
  library(data.table)
  
  print(is.data.frame(TCGA_RNA))
  ####################################
  ####### E X P R E S S I O N ########
  ####################################
  
  # rotate and recreate expression matrix
  print("Transposing expression matrix...")
  
  #TCGA_RNA <- data.table::transpose(TCGA_RNA)
  #colnames(TCGA_RNA) <- rownames(TCGA_RNA)
  #rownames(TCGA_RNA) <- colnames(TCGA_RNA)
  
  colnames(TCGA_RNA)[1] <- "CELLLINE"
  names <-  TCGA_RNA[,1]
  TCGA_RNA.T <- as.data.frame(as.matrix(t(TCGA_RNA[,-1])))
  colnames(TCGA_RNA.T) <- names
  TCGA_RNA <- TCGA_RNA.T
  
  # str_replace(string,"[|]", "..(")
  print("Converting GENE IDs to CCLE format...")
  colnames(TCGA_RNA) <- sapply(stringr::str_replace(colnames(TCGA_RNA),"[|]", "..ENSG000"),function(x) (x[1]))
 
  TCGA_RNA[is.na(TCGA_RNA)] <- 0
  print("Saving final data-sets...")

  TCGA_RNA <- cbind("CELLLINE" = rownames(TCGA_RNA), TCGA_RNA)
  rownames(TCGA_RNA) <- NULL
  
  TCGA_expr  <- TCGA_RNA 
  save(TCGA_expr,file="TCGA_expr.RData")
  #load ("TCGA_expr.RData")
  return()
}