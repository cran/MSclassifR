#############################################
## Variable selection SelectionVar function #
#############################################

SelectionVar <- function(X,
                         Y,
                         MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA"),
                         MethodValidation = c("cv", "repeatedcv", "LOOCV"),
                         PreProcessing = c("center","scale","nzv","corr"),
                         Metric = c("Kappa", "Accuracy"),
                         Sampling = c(NULL, "up", "down", "smote"),
                         NumberCV = NULL,
                         RepeatsCV = NULL,
                         Sizes,
                         Ntree = 1000,
                         threshold = 0.01,
                         ncomp.max = 10){

  
  
  # Methode for preprocess the data (method)
  methodProcess <- PreProcessing

  # Method for training data
  methodValid <- MethodValidation
  
  # Define the method
  SamplingM <- match.arg(Sampling)
  
  switch(SamplingM,
         
  ## Null
  #NULL = {Y <- Y
   #       X <- X}
  
  ## Up
 
  "up" = {
    upTrain <- caret::upSample(x = X,
                             y = Y,
                             list = TRUE)
  
    X <- upTrain$x
    Y <- factor(upTrain$y)
  },
  
  
  ## SMOTE
  "smote" = {
      mozv <- colnames(X)
      
      dataSMOTE <- data.frame(Y,X)
      
      Smoted <- UBL::SmoteClassif(Y~., dataSMOTE, C.perc = "balance")
      
      X <- Smoted[,-1]
      X[is.na(X)] <- 0
    
      colnames(X) = mozv
      Y <- factor(Smoted$Y)
    },
  
  ## Down
  "down" = {
    DownTrain <- caret::downSample(x = X,
                                   y = Y, list = TRUE)
    X <- DownTrain$x
    Y <- factor(DownTrain$y)
  }
  
)


  # Define the method
  method <- match.arg(MethodSelection)
  
  switch(method,

    # RF

    "RFERF" = {

  message("Selection variables with RFE and RF method")
  
  # define the control using a random forest selection function
  control <- caret::rfeControl(functions = caret::rfFuncs,
                               method = methodValid,
                               number = NumberCV,
                               repeats = RepeatsCV)
  
  # run the RFE algorithm
  resultmodel <- caret::rfe(X,
                            Y,
                            rfeControl = control,
                            sizes = Sizes,
                            preProc = PreProcessing,
                            metric = Metric)
  
  #Keeping the peaks giving the best average metric on the folds 
  if (Metric=="Accuracy"){
    nbv=resultmodel$results$Variables[which.max(resultmodel$results$Accuracy[1:(nrow(resultmodel$results)-1)])]
  }
  if (Metric=="Kappa"){
    nbv=resultmodel$results$Variables[which.max(resultmodel$results$Kappa[1:(nrow(resultmodel$results)-1)])]
  }
  
  smoz=unique(resultmodel$variables$var[resultmodel$variables$Variables==nbv])

  results <- (list(result=resultmodel,
              "sel_moz" = sort(smoz)))

  },

    # Logistic regression

    "RFEGlmnet" = {
  
    message("Selection variables with RFE and Glmnet method")
   
  # define the control using a glmnet selection function
  control <- caret::rfeControl(functions = caret::caretFuncs,
                               #lrFuncs
                               method = methodValid,
                               number = NumberCV,
                               repeats = RepeatsCV)
  
  ## Try to fix error task error with specific grid
  #glmnGrid <- expand.grid(alpha = c(0.05, seq(.1, 1, by = 0.5)),
  #                       lambda = c(.001, .01, .1, 1))
  
  
  # run the RFE algorithm
  MyPerf.rfeglmnet <- try(resultmodel <- caret::rfe(X,
                            Y,
                            rfeControl = control,
                            sizes = Sizes,
                            method = "glmnet",
                            preProc = PreProcessing,
                            metric = Metric), silent = T)#,
                            #tuneGrid = glmnGrid)
  
  if(is.character(MyPerf.rfeglmnet)){ 
    warning("glmnet: error probably due to not enough observations by class, the argument MethodValidation is changed by LOOCV.")
    control <- caret::rfeControl(functions = caret::caretFuncs,
                                 #lrFuncs
                                 method = "LOOCV",
                                 number = NumberCV,
                                 repeats = RepeatsCV)
    resultmodel <- caret::rfe(X,
                              Y,
                              rfeControl = control,
                              sizes = Sizes,
                              method = "glmnet",
                              preProc = PreProcessing,
                              metric = Metric)
  }
  
  #Keeping the peaks giving the best average metric on the folds 
  if (Metric=="Accuracy"){
    nbv=resultmodel$results$Variables[which.max(resultmodel$results$Accuracy[1:(nrow(resultmodel$results)-1)])]
  }
  if (Metric=="Kappa"){
    nbv=resultmodel$results$Variables[which.max(resultmodel$results$Kappa[1:(nrow(resultmodel$results)-1)])]
  }
  
  smoz=unique(resultmodel$variables$var[resultmodel$variables$Variables==nbv])
  
  results <- (list(result=resultmodel,
                   "sel_moz" = sort(smoz)))
  },

    #sPLSDA
    "sPLSDA" = {

      message("Selection variables with sPLSDA method")

  ## Select good method for validation
  if (MethodValidation == "LOOCV") {methodValid <- "loo"}
  if (MethodValidation == "LOOCV") {RepeatsCV <- 1}
  if (MethodValidation == "cv") {methodValid <- "Mfold"}
  if (MethodValidation == "cv") {RepeatsCV <- 1}
  if (MethodValidation == "repeatedcv") {methodValid <- "Mfold"}
  if (MethodValidation == "repeatedcv" & is.null(RepeatsCV)) {RepeatsCV <- 1}

  # Select good method for pre process the data
  if (length(grep("scale",PreProcessing)) != 0) {SCALEm <- TRUE} else {SCALEm <- FALSE}
  if (length(grep("nzv",PreProcessing))!= 0) {NZVm <- TRUE} else {NZVm <- FALSE}
  if (length(grep("corr",PreProcessing)) != 0) {warning("The corr preprocessing method is not used with sPLSDA method")}
  if (length(grep("center",PreProcessing)) != 0) {warning("The center preprocessing method is not used with sPLSDA method")}

  if (ncomp.max>(nrow(X)-1)){ncomp.max=nrow(X)-1;}
      
  # Estimate a PLS-DA
  MyResult.plsda <- mixOmics::plsda(X,
                                    Y,
                                    ncomp = ncomp.max,
                                    scale = SCALEm,
                                    near.zero.var = NZVm)

  # Perf function from mixOmics package
  # suggest more for folds and nrepeats
  MyPerf.plsda <- try(mixOmics::perf(MyResult.plsda,
                                 validation = methodValid,
                                 folds = NumberCV,
                                 nrepeat = RepeatsCV,
                                 progressBar = FALSE), silent = T)

  if(is.character(MyPerf.plsda)){ 
    warning("plsa: the system is singular, the argument MethodValidation is replaced by LOOCV.")
    while(is.character(MyPerf.plsda)){
          MyPerf.plsda <- try(mixOmics::perf(MyResult.plsda,
                                   validation = "loo",
                                   folds = NumberCV,
                                   nrepeat = RepeatsCV,
                                   progressBar = FALSE), silent = T)
          ncomp.max=ncomp.max-1;
          MyResult.plsda <- mixOmics::plsda(X,Y,ncomp = ncomp.max,scale = SCALEm,near.zero.var = NZVm)
    }
    MyPerf.plsda <- mixOmics::perf(MyResult.plsda,
                                       validation = "loo",
                                       folds = NumberCV,
                                       nrepeat = RepeatsCV,
                                       progressBar = FALSE)
    warning(paste0(c("plsa: the system is singular, the maximal number of components is replaced by ",ncomp.max,"."),collapse=""))
  }
  

  # Keep optimal number of components for sPLS-DA
  adiff = abs(diff(as.numeric(MyPerf.plsda$error.rate$BER[,1]))) > threshold
  k = 1
  while(isTRUE(adiff[k])){k = k+1;}
  ncomp <- k
  if (ncomp<2){ncomp=2;}

  # tune sPLS-DA to find optimal parameters
  tune.splsda.srbct <- try(mixOmics::tune.splsda(X,Y, ncomp = ncomp,
                                                 validation = methodValid,
                                                 folds = NumberCV,
                                                 dist = "max.dist",
                                                 progressBar = FALSE,
                                                 measure = "BER",
                                                 test.keepX = Sizes), silent = TRUE)
  if(is.character(tune.splsda.srbct)){
    warning("splsa: the system is singular, the argument MethodValidation is changed by LOOCV.")

    tune.splsda.srbct <- mixOmics::tune.splsda(X,Y, ncomp = ncomp,
                                               validation = "loo",
                                               folds = NumberCV,
                                               dist = "max.dist",
                                               progressBar = FALSE,
                                               measure = "BER",
                                               test.keepX = Sizes)
  }

  # Keep optimal number of mass-over-charge values per components
  select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]

  # Train sPLS-DA model with optimal parameters
  splsda.train <- mixOmics::splsda(X,
                                   Y,
                                   scale = SCALEm,
                                   ncomp = ncomp,
                                   keepX = select.keepX,
                                   near.zero.var = NZVm)

  ## Keep best parameters
  comp = ncomp
  data_splsda = splsda.train

  ## Vector of components in the model
  compo <- seq(1:comp)

  ## loop to collect the data
  cp <- list() # Empty list
  variables_keep <- list() # Empty list

  ## loop to extract variable for each component

  for (i in 1:length(compo)){

    ## Keep the data plot
    cp[[i]] <- mixOmics::plotLoadings(data_splsda, comp = compo[[i]],
                                      title = paste("comp",compo[[i]], sep = "_"),
                                      contrib = 'max', method = 'mean', plot=FALSE)

    ## Extract variable in in a vector for each cp
    variables_keep[[i]] <- data.frame(X_var = row.names(cp[[i]]),
                                      Group = cp[[i]]$GroupContrib,
                                      Importance = cp[[i]]$importance)

  }

  ## Rename variable
  names(cp) <- paste("comp_raw", compo, sep = "_")
  names(variables_keep) <- paste("comp", compo, sep = "_")

  ## Remove doublons
  UniqueVariables_keep <- do.call("rbind", variables_keep)
  UniqueVariables_keep <- UniqueVariables_keep[!duplicated(UniqueVariables_keep[,1]),]

  ## Results
  # Final selection of discriminant mass-over-charge values
  results <- list("Raw_data" = cp,
                  "selected_variables" = UniqueVariables_keep,
                  "sel_moz" = sort(UniqueVariables_keep$X_var))

  return(results)

  },

    # VSURF
    
    "VSURF" = {
      
      message("Selection variables with VSURF method")
    resultsmodel <- VSURF::VSURF(x=X, y=Y, ntree = Ntree,
                                 nfor.thres = 50,
                                 nfor.interp = 50,
                                 nfor.pred = 50, nsd=100)

    Vectormoz <- resultsmodel[["varselect.pred"]]
    sel_moz <- colnames(X)[Vectormoz]

    results <- (list(result=resultsmodel,
                     sel_moz=sort(sel_moz)))
  }

  )

  return(results)

}



