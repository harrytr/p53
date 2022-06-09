sigQC_runs <- function()
  
{
  library(sigQC)
  `%notin%` <- Negate(`%in%`)
  setwd("C:/Users/Harry/OneDrive - Nexus365/Desktop/shinyapp")
  wd = getwd()
  print(wd)
  
  regulon_genes <- as.data.frame(read_csv("TP53_regulons.csv"), header  = TRUE)
  regulon_genes <- regulon_genes$x

  
  
  
  setwd(paste0(wd,"//sigQC_data"))
  wd = getwd()
  
  files <- list.files(wd)

  mRNA <- list()
  gene_sigs_list <- list()
  
  tp53_sig <- c("CDC20","CENPA","KIF2C","PLK1") # CDC20,CENPA,KIF2C,PLK1
  
  temp_dir <- ""

  
  for (i in 1: length(files)) {
    
    if (str_sub(files[[i]], start = -4, end = -1) == ".csv") {
      print(files[[i]])
      
      mRNA <- as.data.frame(read_csv(files[[i]]), header  = TRUE)
      print(str_sub(files[[i]], start = 1, end = 4))
      if (str_sub(files[[i]], start = 1, end = 4) != "CCLE") {
        colnames(mRNA)[1] <- "Gene"
        write.csv(mRNA,"mRNA.csv")
        mRNA <- as.data.frame(read_csv("mRNA.csv"), header  = TRUE)
      }


      gene_sigs_list$regulon <- which(mRNA$Gene %in% regulon_genes)
      mRNA[[i]] <- mRNA[,-2]
      names(mRNA[[i]]) <- files[[i]]
    
      temp_dir <- paste0(getwd(),"/",str_sub(files[[i]], start = 1, end = -5))
      print(temp_dir)
      dir.create(path = temp_dir)
      
      setwd( temp_dir)
      
      try(
      g <- make_all_plots(gene_sigs_list, mRNA,  showResults = FALSE, doNegativeControl = FALSE , numResampling = 0, out_dir = temp_dir)
      )
      print(g)
    }
    setwd(wd)
  }
  plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
  plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
  file.copy(from=plots.png.paths, to=temp_dir)
  dev.off()
}
