
MNR_predict <-function(inputs_dir, perm,disease_filename, reg_type, resDir, screening_method, key)  {
  # install.packages(c("glmnet","Matrix","parallel","doParallel","foreach","stats","utils","matrixTests","graphics"))
  #  install.packages(c("testthat","knitr","rmarkdown","plotly","htmltools","kableExtra","DT"))
  #  setwd("C:/Users/Harry/Documents/OneDrive - Nexus365/Desktop/shinyapp/results")
  #Load the library
  library(renoir)
  #Set some directories
  #Input file must be in the directory defined here
  #for the following code to be run without changes
  data = readRDS(paste0(inputs_dir,"//regression.rds"))
  o2 <- sapply(strsplit(colnames(data), split='..', fixed=TRUE),function(x) (x[1]))
  
  
  
  res = readRDS(paste0(inputs_dir,"//best_model.rds"))
  
  # now only select the features from the optimal model in the prediction matrix : 
  optimal_features <- as.data.frame(read.csv(paste0(inputs_dir,"//optimal_features.csv"), header=TRUE))
  optimal_features <- optimal_features$features[2:nrow(optimal_features)]
  
  genes <- paste0(optimal_features,"\\s*?\\.{2}ENSG000",collapse="|")
  genes_e <- paste0("^(", genes, ")")
  data <- data %>%  dplyr::select(matches(genes_e))

  write.csv(data,"sosta.csv")
  print(setdiff(optimal_features,o2))
  analysis.type = reg_type
  #-------------------------------------------------------------------------------#
  #Set X
  
  #Shape data as needed
  dataX = as.matrix(x = data[,-(1),drop=F])
  
  
  #test is numeric
  is.numeric(dataX)
  #attach row names
  rownames(dataX) = data$Variant_Classification
  
  #Check NAs
  #test = sum(colSums(is.na(dataX)))
  
  #-------------------------------------------------------------------------------#
  #Set Y
  
  #Shape data as needed
  
  dataY = as.factor(data[,1])
  table(dataY)
  
  # dataY = droplevels(data[,3,drop=T])
  
  
  #-------------------------------------------------#
  if(identical(analysis.type, "binomial")){
    #BINOMIAL RESPONSE TYPE
    #Response type
    resp.type = "binomial"
    
    #Check which is Missense Mutation
    # index = which(dataY == "Missense_Mutation")
    #Create the Y
    
    dataY = ifelse(stringr::str_detect(dataY,key),1,0)
    #dataY = ifelse(dataY %in% key , yes = 1, no = 0)
    
    #attach row names
    names(dataY) = data$Variant_Classification
    
    #Update x
    dataX = dataX[names(dataY),,drop=F]
    
    
  } else if(identical(analysis.type, "multinomial")){
    #Response type
    resp.type = "multinomial"
    
    #Check data
    #table(dataY)
    
    #Remove In_Frame_Ins (samples numbers too low)
    #
    
    data.sub = data[data$Variant_Classification != "In_Frame_Ins",,drop=F]
    
    
    #Shape data as needed
    
    dataY = droplevels(as.factor(data.sub[,3,drop=T]))
    
    #Check new data
    table(dataY)
    
    #attach row names
    names(dataY) = data.sub$CELLLINE
    
    dataX = dataX[names(dataY),,drop=F]
    
  }
  
  get_typemeasure = function(resp.type){
    type.measure = switch(EXPR = resp.type,
                          # gaussian    = c("mse", "deviance", "mae"),
                          gaussian    = c("mse", "mae"),
                          binomial    = c("mse", "deviance", "class", "auc", "mae"),
                          multinomial = c("mse", "deviance", "class", "mae"));
    return(type.measure);
    
  }
  #as factor
  dataY = as.factor(dataY)
  #-------------------------------------------------------------------------------#
  #Train And Test
  
  #Some parameters to set
  #a) number of repeats: should be set >= 10.
  # Try some numbers. I would use 100, but if too slow you should decrease it.
  n = perm
  #b) number of grid points: should be >=3.
  # 3 is already ok, if execution is fast on your pc you could use 5.
  n.grid.points = 3;
  
  #Write your output directory
  saveRDS(dataY, file = "dataY.rds")
  print("Running renoir to predict on unseen data...")

  
  #Extract the best model
  #a) Select the configuration to use
  best.config = c("min", "1se")[1]
  
  #b) Select the measure type
  type.measure = c("mse", "deviance", "class", "auc", "mae")[3]
  #c) Plot the results
  #plot(object = renoir, best.config = best.config, type.measure = type.measure, nc = length(levels(dataY)))
  

  #d) Select the best model
  best_model = bm.renoir(object = res, best.config = best.config, type.measure = type.measure)
  #e) Probability of recruitment for each feature
  #View(best_model$precruit)
  print("Assessing the performance of the trained best model...")
  accuracy = assess.renoir.trained(
    object = best_model$model,
    dataX,#new prediction matrix 
    dataY, #true data
    best.config = "min"
    #best.config = "1se"#try with this as well, so you have an idea of goodness 
  )
  
  
  saveRDS(accuracy, file = "GLM_testing_results.rds")
 
}
