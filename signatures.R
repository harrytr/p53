signatures <- function(csv) {
  
  library(UpSetR)

  df <- as.data.frame(read.csv(csv, header=TRUE))
  # CCLE_BREAST_DEL/WT	CCLE_ALL_MT/WT	CCLE_BREAST_MT/WT	CCLE_ALL_DEL/WT	TCGA_BRCA_WT_MT	METABRIC_WT/MT
  
  
  CCLE = as.list(unique(as.character(df$CCLE_ALL_MT_WT))) 
  CCLE_BREAST = as.list(unique(as.character(df$CCLE_BREAST_MT_WT)))
  METABRIC = as.list(unique(as.character(df$METABRIC_WT_MT)))
  TCGA_BRCA = as.list(unique(as.character(df$TCGA_BRCA_WT_MT)))
  CCLE_BREAST_DEL = as.list(unique(as.character(df$CCLE_BREAST_DEL_WT)))
  CCLE_ALL_DEL = as.list(unique(as.character(df$CCLE_ALL_DEL_WT)))

  test_bed <- list(CCLE = as.list(unique(as.character(df$CCLE_ALL_MT_WT))), 
                   CCLE_BREAST = as.list(unique(as.character(df$CCLE_BREAST_MT_WT))),
                   METABRIC = as.list(unique(as.character(df$METABRIC_WT_MT))),
                   TCGA_BRCA = as.list(unique(as.character(df$TCGA_BRCA_WT_MT))),
                   CCLE_BREAST_DEL = as.list(unique(as.character(df$CCLE_BREAST_DEL_WT))),
                   CCLE_ALL_DEL = as.list(unique(as.character(df$CCLE_ALL_DEL_WT)))
                    )
  print(upset(fromList(test_bed), order.by = "freq", nsets = 6))


  print(Reduce(intersect,list(METABRIC, TCGA_BRCA, CCLE_BREAST_DEL,CCLE_BREAST,CCLE)))
}