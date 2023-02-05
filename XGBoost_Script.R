    
	   
    
    library(xgboost)
    library(caret)
    library(readxl)
    library(dplyr)
    library(ggplot2)
    library(randomForest)
    library(varImp)
    library(ranger)
    library(janitor)
    library(Boruta)
    library(datasets)
    library(caret)
    library(glmnet)
    library(ggplot2)
    library(data.table)
    print("XGboost training...")


    dataXX <- dataXX %>% dplyr::mutate(Variant_Classification = case_when(str_detect(Variant_Classification,"_Frame_Shift_Del")  ~  1,
                                                       str_detect(Variant_Classification,"_Frame_Shift_Ins")  ~  2,
                                                       str_detect(Variant_Classification,"_In_Frame_Del")  ~  3,
                                                       str_detect(Variant_Classification,"_In_Frame_Ins")  ~  4,
                                                       str_detect(Variant_Classification,"_Missense_Mutation")  ~  5,
                                                       str_detect(Variant_Classification,"_Nonsense_Mutation")  ~  6,
                                                       str_detect(Variant_Classification,"_Splice_Site")  ~  7,
                                                       str_detect(Variant_Classification,"_Fusion")  ~  8,
                                                       str_detect(Variant_Classification,"_Splice_Region")  ~  9,
                                                       str_detect(Variant_Classification,"_Translation_Start_Site")  ~  10))
    
    dataY <- as.factor(dataXX[,1])
    names(dataY) = data$Variant_Classification
    table(dataY)
    cluster <- dataXX[,1]
    clustering_data <- dataXX[,-c(1)]
    clustering_data$cluster <- cluster
    labels <- clustering_data$cluster
    training.samples <- clustering_data$cluster %>%
      createDataPartition(p = 0.7, list = FALSE)
    
    train.data  <- clustering_data[training.samples, ]
    test.data <- clustering_data[-training.samples, ]
    
    print(table(train.data$cluster))
    print(table(test.data$cluster))
    Sys.sleep(4)
    
    train.data <- train.data %>% select(-cluster)
    test.data <- test.data %>% select(-cluster)
    
    labels <- labels[training.samples]
    
    labels <- as.numeric(labels)
    labels <- labels - 1;
    
    
    bstSparse <- xgboost(data = data.matrix(train.data), label=labels, booster = "gbtree",
                         nthread = 8, max_depth = 15,nrounds  = 100,  num_class = 10, objective = "multi:softmax")
    
    print(bstSparse)
    
    importance_matrix <- xgb.importance(model = bstSparse)
    print(importance_matrix)
    write.csv(importance_matrix,"XGBOOST_Selected.csv")
    xgb.plot.importance(importance_matrix = importance_matrix)
    pred <- predict(bstSparse, data.matrix(test.data))
    
    print(pred)
    labels <- as.numeric(clustering_data$cluster)
    labels <- labels[-training.samples]
    labels <- as.numeric(labels)
    labels <- labels - 1;
    print(labels)
    err <- mean(pred != labels)
    print(paste("test-error=", err))