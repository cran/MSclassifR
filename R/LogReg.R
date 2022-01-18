
LogReg <- function(X, moz, Y, number = 2, repeats = 2, kind="linear"){

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

  ## Set up training control
  fit.control <- caret::trainControl(method = "repeatedcv",
                                     number = number,
                                     repeats = repeats,
                                     search="random",
                                     classProbs = T)

  ### Estimation ###
  if (kind=="linear"){
      modelCV <- caret::train(x = DFnnet[, names(DFnnet) != "Y_target"], y = DFnnet$Y_target,
                          method = "multinom", trControl = fit.control, maxit=1000,
                          MaxNWts=84581, preProcess=c("center","scale"), trace = FALSE)
  }
  if (kind=="nnet"){
       modelCV <- caret::train(x = DFnnet[, names(DFnnet) != "Y_target"], y = DFnnet$Y_target,
                           method = "nnet", trControl = fit.control, maxit=1000,
                           MaxNWts=84581, preProcess=c("center","scale"), trace = FALSE)
  }
  if (kind=="rf"){
       modelCV <- caret::train(x = DFnnet[, names(DFnnet) != "Y_target"], y = DFnnet$Y_target,
                             method = "rf", trControl = fit.control, maxit=1000,
                             MaxNWts=84581, preProcess=c("center","scale"), trace = FALSE)
  }
  if (kind=="svm"){
       modelCV <- caret::train(x = DFnnet[, names(DFnnet) != "Y_target"], y = DFnnet$Y_target,
                            method = "svmLinear2", trControl = fit.control, maxit=1000,
                            MaxNWts=84581, preProcess=c("center","scale"), trace = FALSE)
  }
  if (kind=="xgb"){
      xgbGrid <- expand.grid(nrounds = c(1,5,10,12),
                           max_depth = c(1,4,10,16,25),
                           eta = c(.1,.2,.3,.4),
                           gamma = 0,
                           colsample_bytree = c(.5,.7,.8,1),
                           min_child_weight = c(1,5),
                           subsample = c(.8,1))
      modelCV <- caret::train(x = DFnnet[, names(DFnnet) != "Y_target"], y = DFnnet$Y_target,
                            method = "xgbTree", trControl = fit.control,
                            preProcess=c("center","scale"), tuneGrid = xgbGrid, verbose = FALSE)
  }

  ## Arrange data to plot
  variable <- NULL
  value <- NULL
  b_1 <- reshape2::melt(modelCV$resample[,-3])
  b1 <- ggplot2::ggplot(data = b_1 , ggplot2::aes(variable, value, color = factor(variable)))
  b1 = b1 + ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::ylim(0, 1)


  stats_global <- data.frame(mean = apply(modelCV$resample[,-3],2, mean),
                             sd = apply(modelCV$resample[,-3],2, sd))

  ConfusionMATRIX <- confusionMatrix(modelCV)

  return(list(train_mod = modelCV, boxplot = b1, conf_mat = ConfusionMATRIX,
              stats_global = stats_global))
}
