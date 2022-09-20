MNR <-function(regression_data,perm,disease_filename, reg_type, resDir, screening_method, key)  {
  #  install.packages(c("glmnet","Matrix","parallel","doParallel","foreach","stats","utils","matrixTests","graphics"))
  #  install.packages(c("testthat","knitr","rmarkdown","plotly","htmltools","kableExtra","DT"))

  #Load the library
  library(renoir)
  
  #Set some directories
  #Input file must be in the directory defined here
  #for the following code to be run without changes
  
  data = readRDS(regression_data)
  
  colnames(data) <- sapply(strsplit(colnames(data), split='..', fixed=TRUE),function(x) (x[1]))
  
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
  print("Running renoir...")
  res = renoir::run(
    x = dataX,
    y = dataY,
    resp.type = resp.type,
    
    #Pre-processing
    filter = TRUE,
    threshold.na = 0.5,
    threshold.variability = 1,
    na.omit = FALSE,
    
    #Training And Testing
    resampling.method = "multi.random",
    
    n = n,
    min.obs = 10 ,#must be equal to nfolds for glmnet
    n.grid.points = n.grid.points,
    balance = TRUE,
    # balance = FALSE, #you need the bug free version of glmnet to run this option
    do.parallel = TRUE,
    num.cores = 1000,
    restore = FALSE,
    learning.method = "glmnet",
    #Screening
    screening.method = screening_method,
    screening.args = list(sam.args = list(assay.type = "array", logged2 = T),
                          limma.args = list(assay.type = "array", logged2 = T)),
    
    #screening.nvar = seq(50, length(colnames(data)) , round(length(colnames(data))/4)),
    screening.nvar=c(50, 100, 150, 238),
    #Confidence intervals
    confidence = 0.95,
    importance.term = "coefficient",
    random.seed = NULL,
    #Log
    log.path = file.path(resDir, "log.txt"),
    verbose = TRUE,
    
    #Outdir
    outdir = resDir,
    filename = NULL,
    
    #Further arguments
    keep.call = FALSE,
    keep.fit = F,
    relax = FALSE);
  
  #Extract the best model
  #a) Select the configuration to use
  best.config = c("min", "1se")[1]
  # 
  # #b) Select the measure type
  type.measure = c("mse", "deviance", "class", "auc", "mae")[3]
  # #c) Plot the results
  # plot(object = renoir, best.config = best.config, type.measure = type.measure, nc = length(levels(dataY)))
  # 
  # 
  # #d) Select the best model
  #best_model = bm(object = res, best.config = best.config, type.measure = type.measure)
  saveRDS(res,paste0(resDir,"/best_model.rds"))
  #invisible(readline(prompt="Press [enter] to continue"))
  
  # #e) Probability of recruitment for each feature
  # View(best_model$precruit)
  renoir::save_plots_2xpage(res, type.measure = c("mse", "deviance", "class", "mae"),  outdir = resDir, filename = disease_filename)
  
  #YOU NEED TO SET THIS UP WITH YOUR VALUES
  # converison.table must be a dataframe with rows as in dataX columns, or with one column 
  # with the features. The other columns can be annotation (e.g. features is entrez id, another column is hgcn symbol)
  
  #feat.index = the index of the columns of conversion.table that contains the features
  
  renoir::create_report(object = res,
                        outdir = resDir, 
                        filename = paste0('report_multir_',disease_filename, "_", resp.type,'_CCLE'),
                        annotation = NULL, 
                        feat.index = 1,
                        tabset = TRUE)
  
  
}