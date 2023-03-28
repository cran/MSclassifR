
LogReg <- function(X,
                   moz,
                   Y,
                   number = 2,
                   repeats = 2,
                   Metric = c("Kappa", "Accuracy", "F1", "AdjRankIndex", "MatthewsCorrelation"),
                   kind="linear",
                   Sampling = c("no", "up", "down", "smote")){
  message("LogReg function according to the following parameters:")
  moz=as.numeric(moz)
  ### 1. Global matrix of data ###

  ## Extract masses of interest with the global matrix
  X_var <- as.numeric(colnames(X))
  X_var <- data.frame(X_var  = X_var)

  ## Create new matrix with value col
  p_Intmatrix <- cbind(X_var,t(X))

  ### 2. cross-reference the data ###

  "X_var" <- NULL
  Variable.selected <- data.frame(X_var = moz)
  cp_cross <- dplyr::inner_join(Variable.selected, p_Intmatrix, by = "X_var", copy = TRUE)

  ## Col names (features)
  X_cross <- cp_cross[,1] #X_var

  ## TARGET ##
  Y_target <- Y
  Y_target <- factor(Y_target, labels = make.names(levels(Y_target))) # for probs
  
  ## Data
  Data_cross <- cp_cross[,-1] # Clean matrix
  Data_cross <- t(Data_cross) # transpose

  ### 3. Dataframe for multinomial regression ###

  ## renames data
  colnames(Data_cross) <- X_cross
  rownames(Data_cross) <- rownames(X)
  DFnnet <- cbind.data.frame(Y_target = factor(Y_target), Data_cross)
  
  # Metrics non implemented in caret package
  
  ## F1-score function
  
  ### For positive class must be the same that target
  #InterestClass <- make.names(InterestClass)
  
  f1 <- function(data, lev = NULL, model = NULL) {
    f1_val <- f1_score(data$pred,data$obs)
    names(f1_val) <- c("F1")
    f1_val
  }
  
  ## If two classes, F1 will be calculated on the class with the minimum of data
  f1_score <- function(predicted, expected, positive.class = names(which.min(table(Y)))) {
    predicted <- factor(as.character(predicted), levels = unique(as.character(expected)))
    expected  <- as.factor(expected)
    cm = as.matrix(table(expected, predicted))
    
    precision <- diag(cm) / colSums(cm)
    recall <- diag(cm) / rowSums(cm)
    f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    
    #Assuming that F1 is zero when it's not possible compute it
    f1[is.na(f1)] <- 0
    
    #Binary F1 or Multi-class macro-averaged F1
    ifelse(nlevels(expected) == 2, f1[positive.class], mean(f1))
  }
  
  ## Adjusted Rank Index function
  
  adjRankIndex <- function(data, lev = NULL, model = NULL) {
    ARI_val <- mclust::adjustedRandIndex(data$pred,data$obs)
    names(ARI_val) <- c("AdjRankIndex")
    ARI_val
  }

  ##Matthews correlation function
  matthewsCorrelation <- function(data, lev = NULL, model = NULL) {
        Mc_val <- mltools::mcc(data$pred,data$obs)
        names(Mc_val) <- c("MatthewsCorrelation")
        Mc_val
  }
  
  ## Set up training control
  

  ### Modify the smote using an other function ##
  if (isTRUE(Sampling == "smote")) {
    
      dataSMOTE <- DFnnet
      Smoted <- UBL::SmoteClassif(Y_target~., dataSMOTE, C.perc = "balance")
      DFnnet <- Smoted 
    
      Metrica <- match.arg(Metric)
      
      switch(Metrica,
             ## FA
             "F1" = {
               message("smote sampling method and F1 metric selected")
               fit.control <- caret::trainControl(method = "repeatedcv",
                                                  number = number,
                                                  repeats = repeats,
                                                  search="random",
                                                  classProbs = T,
                                                  summaryFunction = f1)},
             
             ## With kappa or accuracy
             "Kappa" = {
               message("smote sampling method and Kappa metric selected")
               fit.control <- caret::trainControl(method = "repeatedcv",
                                                  number = number,
                                                  repeats = repeats,
                                                  search="random",
                                                  classProbs = T)
             },
             
             "Accuracy" = {
               message("smote sampling method and accuracy metric selected")
               fit.control <- caret::trainControl(method = "repeatedcv",
                                                  number = number,
                                                  repeats = repeats,
                                                  search="random",
                                                  classProbs = T)
             },
             
             ### With AdjRankIndex metric
             "AdjRankIndex" = {
               message("smote sampling method and AdjRankIndex metric selected")
               fit.control <- caret::trainControl(method = "repeatedcv",
                                                  number = number,
                                                  repeats = repeats,
                                                  search="random",
                                                  classProbs = T,
                                                  summaryFunction = adjRankIndex)},
             
             ### With Matthews correlation
             "MatthewsCorrelation" = {
               message("smote sampling method and MatthewsCorrelation metric selected")
               fit.control <- caret::trainControl(method = "repeatedcv",
                                                  number = number,
                                                  repeats = repeats,
                                                  search="random",
                                                  classProbs = T,
                                                  summaryFunction = matthewsCorrelation)}
      )
  }
  

  
  else {
    
      Metricb <- match.arg(Metric)
      switch(Metricb,
      ## FA
      "F1" = {
        if(length(Sampling)>1){Sampling="no"}
        if(Sampling=='no'){Sampling=NULL
        message("no sampling method and accuracy metric selected")}else{message(paste(Sampling)," method and F1 metric selected")}
        fit.control <- caret::trainControl(method = "repeatedcv",
                                           number = number,
                                           repeats = repeats,
                                           search="random",
                                           classProbs = T,
                                           summaryFunction = f1,
                                           sampling = Sampling)},
      
      ## With kappa or accuracy
      "Kappa" = {
        if(length(Sampling)>1){Sampling="no"}
        if(Sampling=='no'){Sampling=NULL
        message("no sampling method and Kappa metric selected")}else{message(paste(Sampling)," method and Kappa metric selected")}
        fit.control <- caret::trainControl(method = "repeatedcv",
                                           number = number,
                                           repeats = repeats,
                                           search="random",
                                           classProbs = T,
                                           sampling = Sampling)
      },
      
      "Accuracy" = {
        if(length(Sampling)>1){Sampling="no"}
        if(Sampling=='no'){Sampling=NULL
        message("no sampling method and accuracy metric selected")}else{message(paste(Sampling)," method and accuracy metric selected")}
        fit.control <- caret::trainControl(method = "repeatedcv",
                                           number = number,
                                           repeats = repeats,
                                           search="random",
                                           classProbs = T,
                                           sampling = Sampling)
      },
      
      ### With AdjRankIndex metric
      "AdjRankIndex" = {
        if(length(Sampling)>1){Sampling="no"}
        if(Sampling=='no'){Sampling=NULL
        message("no sampling method and AdjRankIndex metric selected")}else{message(paste(Sampling)," method and AdjRankIndex metric selected")}
        fit.control <- caret::trainControl(method = "repeatedcv",
                                           number = number,
                                           repeats = repeats,
                                           search="random",
                                           classProbs = T,
                                           summaryFunction = adjRankIndex,
                                           sampling = Sampling)},
      
      ### With Matthews correlation
      "MatthewsCorrelation" = {
        if(length(Sampling)>1){Sampling="no"}
        if(Sampling=='no'){Sampling=NULL
        message("no sampling method and MatthewsCorrelation metric selected")}else{message(paste(Sampling)," method and MatthewsCorrelation metric selected")}
        fit.control <- caret::trainControl(method = "repeatedcv",
                                           number = number,
                                           repeats = repeats,
                                           search="random",
                                           classProbs = T,
                                           summaryFunction = matthewsCorrelation,
                                           sampling = Sampling)}
      )
    
  }
  
 
  
  #####
  
  
  ### Estimation ###
  if (kind=="linear"){
    message("estimation with linear method")
      modelCV <- caret::train(Y_target~.,
                              data = DFnnet,
                              method = "multinom",
                              trControl = fit.control,
                              maxit=1000,
                              MaxNWts=84581,
                              preProcess=c("center","scale"),
                              trace = FALSE,
                              metric = Metric)
  }

  row.names(DFnnet) = NULL
  if (kind=="nnet"){
    message("estimation with nnet method")
       modelCV <- caret::train(Y_target~.,
                               data = DFnnet,
                               method = "nnet",
                               trControl = fit.control,
                               maxit=1000,
                               MaxNWts=84581,
                               tuneLength = 10,
                               preProcess=c("center","scale"),
                               trace = FALSE,
                               metric = Metric)
  }
  if (kind=="rf"){
    message("estimation with rf method")
       modelCV <- caret::train(Y_target~.,
                               data = DFnnet,
                               method = "rf",
                               trControl = fit.control,
                               maxit=1000,
                               MaxNWts=84581,
                               tuneLength = 10,
                               preProcess=c("center","scale"),
                               trace = FALSE,
                               metric = Metric)
  }
  if (kind=="svm"){
    message("estimation with svm method")
       modelCV <- caret::train(Y_target~.,
                               data = DFnnet,
                               method = "svmLinear2",
                               trControl = fit.control,
                               maxit=1000,
                               MaxNWts=84581,
                               tuneLength = 10,
                               preProcess=c("center","scale"),
                               trace = FALSE,
                               metric = Metric)
  }
  if (kind=="xgb"){
    message("estimation with xgbTree method")
      ## Grid for xgbTree
      xgbGrid <- expand.grid(nrounds = c(1,5,10,12),
                           max_depth = c(1,4,10,16,25),
                           eta = c(.1,.2,.3,.4),
                           gamma = 0,
                           colsample_bytree = c(.5,.7,.8,1),
                           min_child_weight = c(1,5),
                           subsample = c(.8,1))
      
      
      modelCV <- caret::train(Y_target~.,
                              data = DFnnet,
                              method = "xgbTree",
                              trControl = fit.control,
                              preProcess=c("center","scale"),
                              tuneGrid = xgbGrid,
                              verbose = FALSE,
                              metric = Metric,
                              verbosity = 0)
    
  }
  
  
  ## Arrange data to plot
  variable <- NULL
  value <- NULL
  b_1 <- suppressMessages(reshape2::melt(modelCV$resample[,-3]))
  b1 <- ggplot2::ggplot(data = b_1 , ggplot2::aes(variable, value, color = variable))
  b1 = b1 + ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::ylim(0, 1)

  Mean_Metric <- dplyr::summarize(dplyr::group_by(b_1, variable), Mean = mean(value, na.rm=TRUE))
  Sd_metric <- dplyr::summarize(dplyr::group_by(b_1, variable), Sd = sd(value, na.rm=TRUE))
  statsGlobal <- merge(Mean_Metric, Sd_metric)
  colnames(statsGlobal) <- c("Metric", "Mean", "Sd")
  
  Confusion.Matrix <- confusionMatrix(modelCV, "none")
  
  return(list(train_mod = modelCV,
              boxplot = b1,
              Confusion.Matrix = Confusion.Matrix,
              stats_global = statsGlobal))
}
