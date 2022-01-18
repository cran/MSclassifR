#############################################
## Variable selection SelectionVar function #
#############################################

SelectionVar <- function(X,
                         Y,
                         MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA"),
                         MethodValidation = c("cv", "repeatedcv", "LOOCV"),
                         PreProcessing = c("center","scale","nzv","corr"),
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
  method <- match.arg(MethodSelection)

  switch(method,

    # RF

    "RFERF" = {

  print("Selection variables with RFE and RF method")
  
  # define the control using a random forest selection function
  control <- caret::rfeControl(functions = caret::rfFuncs,
                               method= methodValid,
                               number = NumberCV,
                               repeats = RepeatsCV)

  # run the RFE algorithm
  resultmodel <- caret::rfe(X,
                        Y,
                        rfeControl= control,
                        sizes = Sizes,
                        preProc = PreProcessing)


  results <- (list(result=resultmodel,
              "sel_moz" = as.numeric(resultmodel[["optVariables"]])))

  },

    # Logistic regression

    "RFEGlmnet" = {
  
  print("Selection variables with RFE and Glmnet method")
   
  # define the control using a glmnet selection function
  control <- caret::rfeControl(functions = caret::caretFuncs,
                               method = methodValid,
                               number = NumberCV,
                               repeats = RepeatsCV)

  # run the RFE algorithm
  resultmodel <- caret::rfe(X,
                            Y,
                            rfeControl = control,
                            sizes = Sizes,
                            method = "glmnet",
                            preProc = PreProcessing)
                           


  results <- (list(result=resultmodel,
                   "sel_moz" = as.numeric(resultmodel[["optVariables"]])))
  },

    #sPLSDA
    "sPLSDA" = {

  print("Selection variables with sPLSDA method")

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
  if (length(grep("corr",PreProcessing)) != 0) {warning("this preprocessing method is not avaible with sPLSDA method")}
  if (length(grep("center",PreProcessing)) != 0) {warning("this preprocessing method is not avaible with sPLSDA method")}

  # Estimate a PLS-DA
  MyResult.plsda <- mixOmics::plsda(X,
                                    Y,
                                    ncomp = ncomp.max,
                                    scale = SCALEm,
                                    near.zero.var = NZVm)

  # Perf function from mixOmics package
  # suggest more for folds and nrepeats
  MyPerf.plsda <- mixOmics::perf(MyResult.plsda,
                                 validation = methodValid,
                                 folds = NumberCV,
                                 nrepeat = RepeatsCV,
                                 progressBar = FALSE)

  # Keep optimal number of components for sPLS-DA
  adiff = abs(diff(as.numeric(MyPerf.plsda$error.rate$BER[,1]))) > threshold
  k = 1
  while(isTRUE(adiff[k])){k = k+1;}
  ncomp <- k

  # tune sPLS-DA to find optimal parameters
  tune.splsda.srbct <- mixOmics::tune.splsda(X,Y, ncomp = ncomp,
                                             validation = methodValid,
                                             folds = NumberCV,
                                             dist = 'max.dist',
                                             progressBar = FALSE,
                                             measure = "BER",
                                             test.keepX = Sizes)

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
                                      contrib = 'max', method = 'mean')

    ## Extract variable in in a vector for each cp
    variables_keep[[i]] <- data.frame(X_var = as.numeric(row.names(cp[[i]])),
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
                  "sel_moz" = as.numeric(UniqueVariables_keep$X_var))

  return(results)

  },

    # VSURF
    
    "VSURF" = {
      
    print("Selection variables with VSURF method")
    resultsmodel <- VSURF::VSURF(X, Y, ntree = Ntree,
                                 nfor.thres = 50,
                                 nfor.interp = 50,
                                 nfor.pred = 50)

    Vectormoz <- resultsmodel[["varselect.pred"]]
    sel_moz <- colnames(X)[Vectormoz]

    results <- (list(result=resultsmodel,
                     sel_moz=sel_moz))
  }


  )

  return(results)

}



