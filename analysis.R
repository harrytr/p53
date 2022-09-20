
MNR <-function(regression_data,perm,disease_filename, reg_type, resDir, screening_method, key)  {
  # install.packages(c("glmnet","Matrix","parallel","doParallel","foreach","stats","utils","matrixTests","graphics"))
  #  install.packages(c("testthat","knitr","rmarkdown","plotly","htmltools","kableExtra","DT"))

  library(renoir)

  
  get_hp = function(id, y){
    
    #Generalised Linear Models with Penalisation
    lambda = 10^seq(3, -2, length=100)
    # alpha = seq(0.1, 0.9, length = 9)
    alpha = seq(0.1, 0.9, length = 5)
    gamma = c(0, 0.25, 0.5, 0.75, 1)
    
    #Random Forest
    ntree = c(10, 50, 100, 250, 500)
    
    #Generalised Boosted Regression Modelling
    eta = c(0.3, 0.1, 0.01, 0.001)
    
    #Support Vector Machines
    cost      = 2^seq(from = -5, to = 15, length.out = 5)
    svm.gamma = 2^seq(from = -15, to = 3, length.out = 4)
    degree    = seq(from = 1, to = 3, length.out = 3)
    #Note that for classification nu must be
    #nu * length(y)/2 <= min(table(y)). So instead of
    #fixing it as
    # nu        = seq(from = 0.1, to = 0.6, length.out = 5)
    #we do
    nu.to = floor((min(table(y)) * 2/length(y)) * 10) / 10
    nu = seq(from = 0.1, to = nu.to, length.out = 5)
    
    #kNN
    k = seq(from = 1, to = 9, length.out = 5)
    
    #Nearest Shrunken Centroid
    threshold = seq(0, 2, length.out = 30)
    
    #hyperparameters
    out = switch(
      id,
      'lasso'              = list(lambda = lambda),
      'ridge'              = list(lambda = lambda),
      'elasticnet'         = list(lambda = lambda, alpha = alpha),
      'relaxed_lasso'      = list(lambda = lambda, gamma = gamma),
      'relaxed_ridge'      = list(lambda = lambda, gamma = gamma),
      'relaxed_elasticnet' = list(lambda = lambda, gamma = gamma, alpha = alpha),
      'randomForest'       = list(ntree = ntree),
      'gbm'                = list(eta = eta, ntree = ntree),
      'linear_SVM'         = list(cost = cost),
      'polynomial_SVM'     = list(cost = cost, gamma = svm.gamma, degree = degree),
      'radial_SVM'         = list(cost = cost, gamma = svm.gamma),
      'sigmoid_SVM'        = list(cost = cost, gamma = svm.gamma),
      'linear_NuSVM'       = list(nu = nu),
      'polynomial_NuSVM'   = list(nu = nu, gamma = svm.gamma, degree = degree),
      'radial_NuSVM'       = list(nu = nu, gamma = svm.gamma),
      'sigmoid_NuSVM'      = list(nu = nu, gamma = svm.gamma),
      'gknn'               = list(k = k),
      'nsc'                = list(threshold = threshold)
    )
    
    return(out)
  }
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
  
  print(is.numeric(dataX))

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
    print("Binomial classification starting...")
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
  # NEW SET UP
  #-------------------------------------------------------------------------------#
  #learning method
  learning.method.id = "elasticnet"

  #sampling for tuning
  sampling.method.id.tuning = "cv"

  #sampling for evaluation
  sampling.method.id.evaluation = "random"

  #metric for tuning
  performance.metric.id.tuning = "acc" #MUST BE ONE METRIC

  #metrics for evaluation
  performance.metric.ids.evaluation = c("acc", "precision") #YOU CAN WRITE HERE MORE THAN ONE METRIC

  #-------------------------------------------------------------------------------#
  #Filter
  filterL = FilterList(
    Filter(id = "na"),
    Filter(id = "intensity"),
    Filter(id = "variability")
  )

  #-------------------------------------------------------------------------------#
  #Screeners
  screeners = ScreenerList(
    Screener(id = "ebayes", parameters = list(assay.type = "array", logged2 = T)),
    Screener(id = "permutation", parameters = list(assay.type = "array", logged2 = T))
  )

  #-------------------------------------------------------------------------------#
  #tuner

  tuner = Tuner(
    id = "grid.search",
    sampler = Sampler(
      method = sampling.method.id.tuning,
      k = 10L,
      n = integer(),
      strata = dataY
    ),
    screener = screeners,
    looper   = Looper(cores = 1L),
    logger   = Logger(verbose = T, level = "INFO")
  )


  #-------------------------------------------------------------------------------#
  #Learner
  learner = Learner(
    tuner      = tuner,
    trainer    = Trainer(id = learning.method.id),
    forecaster = Forecaster(id = learning.method.id),
    scorer     = ScorerList(Scorer(id = performance.metric.id.tuning)),
    selector   = Selector(id = learning.method.id),
    recorder   = Recorder(id = learning.method.id),
    marker     = Marker(id = learning.method.id),
    logger     = Logger(level = "ALL")
  )


  #-------------------------------------------------------------------------------#
  #ScorerList
  sl = list()
  for(i in performance.metric.ids.evaluation){
    sl[[i]] = Scorer(id = i)
  }
  sl = ScorerList(sl)

  #Evaluator

    
    evaluator = Evaluator(
      #Sampling strategy: stratified random sampling without replacement
      sampler = Sampler(               
        method = "random",             
        k = 10L,                       
        strata = dataY,                    
        N = as.integer(length(dataY))      
      ),
      
    #Performance metric
    #Each chosen metric must be in a different Scorer, so if you select 3 metrics you will have to add another Scorer in the ScorerList
    scorer  = sl,

    #looper  = Looper(cores = 1L)#sequential Loop
    looper  = Looper(cores = 1000L)#parallel Loop
  )

  #-------------------------------------------------------------------------------#
  # Train And Test
  #-------------------------------------------------------------------------------#

  saveRDS(dataY, file = "dataY.rds")
  print("Running renoir...")
  res <- renoir(
    filter = filterL,

    #Number of considered training set sizes
    npoints = 3,
    # ngrid,
    ##Minimum number of samples in the grid
    nmin = round(nrow(dataX)/2),
    #nmin = 5,

    #Loop
    looper = Looper(),

    #Store
    filename = "renoir",
    outdir   = resDir,
    restore  = TRUE,

    #Learn
    learner   = learner,

    #Evaluate
    evaluator = evaluator,

    #Log
    logger    = Logger(path = file.path(resDir, "log.txt"), level = "ALL", verbose = T),

    #Data for training
    hyperparameters = get_hp(id = learning.method.id, y = y),
    x         = dataX,
    y         = dataY,
    weights   = NULL,
    offset    = NULL,
    resp.type = resp.type,

    #Free space
    rm.call = FALSE,
    rm.fit  = FALSE,

    #Group results
    grouping = TRUE,

    #Screening
    screening = NULL,

    #Remove call from trainer to reduce space
    keep.call = F
  )

  saveRDS(res,paste0(resDir,"/best_model.rds"))
  
  renoir:::create_report(object = res[[1]],
                        outdir = resDir,
                        filename = paste0('report_multir_',disease_filename, "_", resp.type,'_CCLE'),
                        report.type = "full",
                        annotation = NULL,
                        feat.index = 1,
                        tabset = TRUE)

  # renoir:::plot.RenoirSummaryTable(
  #   x = res[res$config == "opt",,drop=F],
  #   measure     = "precision", 
  #   set         = "full", 
  #   interactive = T, 
  #   add.boxplot = F,
  #   add.scores  = F,
  #   add.best    = F,
  #   key         = c("id", "config")
  # )
  # 

}
