#############################################
## Variable selection SelectionVar function #
#############################################

SelectionVar <- function(X,
                         Y,
                         MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA", "mda", "cvp"), 
                         MethodValidation = c("cv", "repeatedcv", "LOOCV"),
                         PreProcessing = c("center", "scale", "nzv", "corr"), 
                         Metric = c("Kappa", "Accuracy"),
                         Sampling = c("no", "up","down", "smote"),
                         NumberCV = NULL, RepeatsCV = NULL, 
                         Sizes,
                         Ntree = 1000,
                         ncores = 2,
                         threshold = 0.01,
                         ncomp.max = 10,
                         nbf=0){
  
  # Methode for preprocess the data (method)
  methodProcess <- PreProcessing
  
  # Method for training data
  methodValid <- MethodValidation
  
  # Define the method
  SamplingM <- match.arg(Sampling)
  
  switch(SamplingM,
         
         "no" = {
           message("No sampling method selected")
           X=X; Y=Y;},
         
         ## Null
         #NULL = {Y <- Y
         #       X <- X}
         
         ## Up
         
         "up" = {
           message("Up sampling method selected")
           upTrain <- caret::upSample(x = X,
                                      y = Y,
                                      list = TRUE)
           
           X <- upTrain$x
           Y <- factor(upTrain$y)
         },
         
         
         ## SMOTE
         "smote" = {
           message("Smote sampling method selected")
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
           message("Down sampling method selected")
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
         },
         
         # mda
         
         "mda" = {
           message("Selection variables with mda method")
           if (nbf>0){
             #creating false peaks to get enough variables with negative importances
             X0=matrix(stats::runif(n = nbf*nrow(X),min = min(X), max = max(X)),nrow = nrow(X),ncol=nbf)
             colnames(X0)=paste("false_",1:ncol(X0),sep="")
             #Adding them to X
             Xn=cbind(X,X0);
           }else{Xn=X;}
           #Estimating random forests
           rf = randomForest::randomForest(x = Xn, Y, ntree = Ntree, 
                                           mtry = floor(sqrt(ncol(X))), importance = TRUE, keep.forest = F)
           vi = randomForest::importance(rf, type = 1, scale = FALSE)
           
           if (nbf>0){
             #Deleting false with positive importances
             vi1 = c(vi[1:ncol(X)],
                     vi[(ncol(X)+1):length(vi)][which(vi[(ncol(X)+1):length(vi)]<0)])
             names(vi1)=c(rownames(vi)[1:ncol(X)],rownames(vi[(ncol(X)+1):length(vi)][which(vi[(ncol(X)+1):length(vi)]<0)]))
           }else{vi1 = c(vi[1:ncol(X)]);names(vi1)=rownames(vi)[1:ncol(X)];}
           Fall = stats::ecdf(vi1)
           imp_neg = vi1[which(vi1 < 0)]
           #Creating distribution under the null by the opposite of negative importances
           imp_null = c(imp_neg, -imp_neg)
           #Estimating pi0 for different quantile values of the null
           q_ext = seq(0.75, 1, by = 0.01)
           pi0 = NULL
           for (i in q_ext) {
             qin = stats::quantile(imp_null, i)
             pi0 = c(pi0, min(Fall(qin)/i, 1))
           }
           if (nbf>0){
             Nfn=sum(vi1[(ncol(X)+1):length(vi1)]<0)
             #Reajusting the estimation on the set of true peaks
             pi0f = (min(pi0)*(ncol(X)+Nfn)-Nfn)/ncol(X)
           }else{pi0f = min(pi0);}
           #Results
           nb_to_sel = floor(ncol(X) * (1 - pi0f))
           sel_moz = names(vi1[1:ncol(X)])[sort(-vi1[1:ncol(X)], index.return = TRUE)$ix][1:nb_to_sel]
           imp_sel = vi1[which(names(vi1) %in% sel_moz)]
           results <- list("nb_to_sel" = nb_to_sel,
                           "sel_moz" = sel_moz, 
                           "imp_sel" = imp_sel)
         },
         
         "cvp" = {
           message("Selection variables with cvp method")
           if (is.null(NumberCV)){NumberCV=2;}
           
           if (nbf>0){
             #creating false peaks to get enough variables with negative importances
             X0=matrix(runif(n = nbf*nrow(X),min = min(X), max = max(X)),nrow = nrow(X),ncol=nbf)
             colnames(X0)=paste("false_",1:ncol(X0),sep="")
             #Adding them to X
             Xn=cbind(X,X0)
           }else{Xn=X;}
           
           cv_vi = vita::CVPVI(Xn, as.numeric(Y), k = NumberCV, ntree = Ntree, 
                               ncores = ncores)
           vi = cv_vi$cv_varim
           if (nbf>0){
             #Deleting false with positive importances
             vi1 = c(vi[1:ncol(X)],
                     vi[(ncol(X)+1):length(vi)][which(vi[(ncol(X)+1):length(vi)]<0)])
             names(vi1)=c(rownames(vi)[1:ncol(X)],rownames(vi[(ncol(X)+1):length(vi)][which(vi[(ncol(X)+1):length(vi)]<0)]))
           }else{vi1 = c(vi[1:ncol(X)]);names(vi1)=rownames(vi)[1:ncol(X)];}
           Fall = stats::ecdf(vi1)
           imp_neg = vi1[which(vi1 < 0)]
           #Creating distribution under the null by the opposite of negative importances
           imp_null = c(imp_neg, -imp_neg)
           q_ext = seq(0.75, 1, by = 0.01)
           pi0 = NULL
           for (i in q_ext) {
             qin = stats::quantile(imp_null, i)
             pi0 = c(pi0, min(Fall(qin)/i, 1))
           }
           if (nbf>0){
             Nfn=sum(vi1[(ncol(X)+1):length(vi1)]<0)
             pi0f = (min(pi0)*(ncol(X)+Nfn)-Nfn)/ncol(X)
           }else{pi0f = min(pi0);}
           nb_to_sel = floor(ncol(X) * (1 - pi0f))
           sel_moz = names(vi1[1:ncol(X)])[sort(-vi1[1:ncol(X)], index.return = TRUE)$ix][1:nb_to_sel]
           imp_sel = vi1[which(names(vi1) %in% sel_moz)]
           results <- list("nb_to_sel" = nb_to_sel,
                           "sel_moz" = sel_moz, 
                           "imp_sel" = imp_sel)
         }
  )
  return(results)
}
