########################
### Predict Category ###
########################

PredictLogReg <-  function (peaks,
                            model,
                            moz,
                            tolerance = 6,
                            toleranceStep = 2,
                            normalizeFun = TRUE,
                            noMatch = 0,
                            Reference = NULL){
  
  if(length(unique(names(peaks))) != length(peaks) && is.null(names(peaks)) == FALSE){stop("Each element of X must have a unique name")}
 
  #Create a subfunction
  predi=function(Peaks,
                 Modell,
                 Moz,
                 Tolerance,
                 NormalizeFun,
                 noMatch){
    
    ## Empty list
    Peak <- list()
    DF.match <- list()
    DF.Peaks <- list()
    prediction <- list()
    diff_mass <- list()
    DF.matchSplit <- list()
    name=rep(NA,length(Peaks))
    method=name
    Tolerance.step <- list()
    
    for (i in 1:length(Peaks)){
      
      if(is.null(Peaks[[i]]@metaData$fullName)){name[i]=Peaks[[i]]@metaData$file;}
      else{name[i]=Peaks[[i]]@metaData$fullName;}
      
      
      # Print Metadata
      print(name[i])
      
      method[i]= Modell$method;
      
      Peak[[i]] <- data.frame(X_var = Peaks[[i]]@mass, # mass over charges to match
                              intensity = Peaks[[i]]@intensity) # intensities
      
      ## Sort the selected mass over charge
      s_moz <- data.frame ("X_var" = unique(sort(as.numeric(Moz))))
      
      ## Match intensities of the peaks to the selected mass over charge values using the Tolerance argument
      DF.match[[i]] <- fuzzyjoin::distance_left_join(s_moz,
                                                     Peak[[i]],
                                                     by = "X_var",
                                                     max_dist = Tolerance)
      
      ### If no match
      if(sum(is.na(DF.match[[i]]$X_var.y)) == nrow(DF.match[[i]])){
        warning("No m/z from peaks object matched with the m/z in the moz objet,
                the tolerance has been increased by steps indicated in toleranceStep argument Da")
        Tolerance.step <- Tolerance
        
        
        while(sum(is.na(DF.match[[i]]$X_var.y)) == nrow(DF.match[[i]])){
          Tolerance.step <- Tolerance.step + toleranceStep  #Tolerance step
          
          DF.match[[i]] <- fuzzyjoin::distance_left_join(s_moz,
                                                         Peak[[i]],
                                                         by = "X_var",
                                                         max_dist = Tolerance.step)
          Tolerance.step <- Tolerance.step
          
        }
        
        message(paste(c("tolerance found for match =", Tolerance.step)))
      }
      
      #diff_mass[[i]] <- abs(DF.match[[i]][,2] - DF.match[[i]][,1])
      #DF.match[[i]] <- cbind.data.frame(DF.match[[i]], diff_mass = diff_mass[[i]])
      #DF.matchSplit[[i]] <- split.data.frame(DF.match[[i]], DF.match[[i]]$X_var.y)
      #DF.matchSplit[[i]] <- do.call("rbind",lapply(DF.matchSplit[[i]], function(x) x[which.min(x[,4]),]))
      #DF.match[[i]] <- dplyr::full_join(DF.matchSplit[[i]], data.frame(X_var.x = as.numeric(moz)), by = "X_var.x")
      #DF.match[[i]] <- DF.match[[i]][order(DF.match[[i]][,1]),][,-4]
      
      ## When no match, intensity values are replaced by the noMatch argument
      DF.match[[i]]$intensity[is.na(DF.match[[i]]$intensity)] <- noMatch
      
      #if several matchs
      td=table(DF.match[[i]]$X_var.x)
      dbl=names(td)[td>1]
      if (length(dbl)>0){
        keepDF=NULL
        for (l in 1:length(dbl)){
          diff_m=abs(DF.match[[i]][DF.match[[i]]$X_var.x==dbl[l],2]-DF.match[[i]][DF.match[[i]]$X_var.x==dbl[l],1])
          keepDF=rbind(keepDF,DF.match[[i]][DF.match[[i]]$X_var.x==dbl[l],][which.min(diff_m),])
        }
        DF.match[[i]]=DF.match[[i]][which(!DF.match[[i]]$X_var.x%in%dbl),]
        DF.match[[i]]=rbind(DF.match[[i]],keepDF)
        DF.match[[i]]=DF.match[[i]][order(DF.match[[i]][,1]),]
      }
      
      DF.Peaks[[i]] <- data.frame(DF.match[[i]]$intensity)
      DF.Peaks[[i]] <- t(DF.Peaks[[i]])
      colnames(DF.Peaks[[i]]) <- as.character(DF.match[[i]]$X_var.x)
      
      # Normalize peak with the maximum intensity value
      if (NormalizeFun){
        norma<-function(x) x/(max(x))
        DF.Peaks[[i]] <- apply(DF.Peaks[[i]],1,norma)
        DF.Peaks[[i]] <- t(DF.Peaks[[i]])
      }
      
      prediction[[i]] <- stats::predict(Modell,DF.Peaks[[i]], type = "prob");
      
      
    }
    
    ## Results
    Results <- do.call("rbind", prediction)
    Results=cbind(name,method,Results)
    rownames(Results)=paste(1:nrow(Results),".",sep="")
    NameStrain <- rep(names(peaks), nlevels(factor(method)))
    if (is.null(names(peaks))){Results$name = rep(c(1:length(peaks)), nlevels(factor(method)))} else {Results$name = rep(NameStrain, nlevels(factor(method)))}
    return(Results)
  }
  
  ## Create a data frame for mass spectra Peaks
  if (!is.null(model$method)){
    model=list(model)
  }
  res=NULL
  for (i in 1:length(model)){
    res1=predi(Peaks=peaks,Modell=model[[i]],Moz= as.numeric(moz),Tolerance=tolerance,
               NormalizeFun=normalizeFun,noMatch=noMatch)
    res=rbind(res,res1)
  }
  ##Merging results using Fisher method
  if (length(model)>1){
    lr=levels(as.factor(res$name))
    if (length(lr)>1){
      for (i in 1:length(lr)){
        #Fisher's method
        res_1=res[res$name==lr[i],]
        suppressWarnings({
          merge_p= apply(res_1[,3:ncol(res_1)],2,function(x)return(metap::sumlog(as.numeric(x))$p))
        })
        #Max vote
        res_1=res[res$name==lr[i],3:ncol(res_1)]
        Max1 = apply(data.frame(res_1), 1, function(x) replace(as.numeric(x), which.max(x), 1))
        Max2 = apply(Max1, 2, function(x) replace(x, x<1, 0))
        VoteTot = ncol(Max2)
        Votedecompte = apply(Max2, 1, sum)
        Votedecompte2 = Votedecompte/VoteTot
        mat_vote = Votedecompte2
        #####"
        #votes=table(unlist(apply(res_1,1,function(x){return(which(as.numeric(x)==max(as.numeric(x))))})))
        #vote_max=which(as.numeric(votes)==max(as.numeric(votes)))
        #mat_vote=rep(0,ncol(res_1))
        #mat_vote[which(as.numeric(votes)==max(as.numeric(votes)))]=1
        #mat_vote=mat_vote/sum(mat_vote)
        
        res=rbind(res,c(lr[i],"comb_fisher",merge_p))
        res=rbind(res,c(lr[i],"max_vote",mat_vote))
      }
    }
  }
  
  #Add a column indicating the predicted category
  pred_cat=NULL
  for (i in 1:nrow(res)){
    pred_cat=c(pred_cat,colnames(res)[3:ncol(res)][which.max(res[i,3:ncol(res)])])
  }
  
  res=cbind(res,pred_cat)
  colnames(res)[ncol(res)]="pred_max_p"
  
  ## Statistcs for analysis
  ## For models
  if (is.null(Reference) == FALSE && length(unique(as.character(Reference))) > 1){
    df_methods <- split.data.frame(res, res$method)
    suppressWarnings({
      RefT = factor(make.names(Reference))
      
      Model_prediction <- lapply(df_methods, function(x) factor((x[,"pred_max_p"]), levels = levels(RefT)))
      
      #Model_prediction <- Model_prediction[names(Model_prediction) %in% c("max_vote", "comb_fisher") == FALSE]
      
      Confusion.Matrix <- lapply(Model_prediction, function(x) try(caret::confusionMatrix(x, RefT), silent = TRUE))
      Confusion.MatrixT <- lapply(Confusion.Matrix, function(x) if(is.character(x) == FALSE){x<-x})
      Confusion.MatrixT <- Confusion.MatrixT[!sapply(Confusion.MatrixT,is.null)]
      
      Confusion.Matrix = Confusion.MatrixT 
    })
    ## More adapted Warning
    #n.levels <- as.numeric(sapply(Model_prediction, nlevels))
    #n.levelsRef <- as.numeric(nlevels(RefT))
    #cond.num <- sum(as.numeric(n.levels%in%n.levelsRef))
    #if(cond.num <=1){message("Be careful, at least one model is not able to predict certain classes")}
    
    Df.classif <- lapply(df_methods, function(x) try(cbind.data.frame(Reference = make.names(Reference), x), silent = TRUE))
    Df.classif <- Df.classif[!sapply(Df.classif,is.character)]
    DiscordTable <- lapply(Df.classif, function(x) x[which(x[,"Reference"]!= x[,"pred_max_p"]),])
    
    #table for discordances strains
    #Discordances.strain <- lapply(DiscordTable, function(x) table(x$name))
    #Discordances.strains <- reshape2::melt(Discordances.strain, id.vars = NULL)
    
    ConfMat.df <- lapply(Confusion.Matrix, function(x) data.frame(x[["table"]]))
    Discord2 <- lapply(ConfMat.df, function(x) x[which(x[,1]!= x[,2]),])
    Incor.Classification <- lapply(Discord2, function(x) x[which(x[,3]>=1),])
    
    # Incorrect classification with frequencies
    suppressMessages(Incorrect.ClassificationFreq <- try(reshape2::melt(Incor.Classification)[,-3], silent = TRUE))
    
    if(is.character(Incorrect.ClassificationFreq)){Incorrect.ClassificationFreq <- Incorrect.ClassificationFreq[!sapply(Incorrect.ClassificationFreq,is.character)]}
    if(is.character(Incorrect.ClassificationFreq)){Incorrect.ClassificationFreq <- t(data.frame(c(0,0,0,0)))}
    colnames(Incorrect.ClassificationFreq) <- c("Prediction", "Reference", "Frequence", "Model")
    
    ConcordTable <- lapply(Df.classif, function(x) x[which(x[,"Reference"]== x[,"pred_max_p"]),])
    
    #table for concordances strains
    #Concordances.strain <- lapply(DiscordTable, function(x) table(x$name))
    #Concordances.strains <- reshape2::melt(Concordances.strain, d.vars = NULL)
    
    Cor.Classification <- lapply(ConfMat.df, function(x) x[which(x[,1] == x[,2]),])
    Cor.Classification <- lapply(Cor.Classification, function(x) x[which(x[,3]>0),])
    
    # Correct classification with frequencies
    suppressMessages(Correct.ClassificationFreq <- reshape2::melt(Cor.Classification)[,-3])
    colnames(Correct.ClassificationFreq) <- c("Prediction", "Reference", "Frequence", "Model")
    
    # Statistics
    Global.stat <- lapply(Confusion.Matrix, function(x) (x[["overall"]]))
    Global.stats <- do.call("data.frame", Global.stat)
    Statistic.param = row.names(Global.stats)
    
    Global.stats <- suppressMessages(reshape2::melt(Global.stat, id.vars = NULL))
    Global.stats <- data.frame(Global.stats, Statistic.param)
    
    Details.stat <- try(lapply(Confusion.Matrix, function(x) cbind.data.frame(x[["byClass"]], "Class" = row.names(x[["byClass"]]))), silent = TRUE)
    suppressMessages(Details.stats <- try(reshape2::melt(Details.stat), silent = TRUE)) 
    
    ## If only two classes
    if(is.character(Details.stat)){
      Details.stat2 <-lapply(Confusion.Matrix, function(x) cbind.data.frame(x[["byClass"]], "Class" = x[["positive"]]))
      suppressMessages(Details.stats2 <- cbind.data.frame(reshape2::melt(Details.stat2)))
      colnames(Details.stats2) <- c("Class", "Statistic.parameter", "Value", "Model")
      Details.stats2$Statistic.parameter <- row.names(data.frame(Details.stat2[1]))
      
      Details.stats = Details.stats2
    }
    # Other statistics
    ## Matthew Correlation for mulitclass
    MatthewCorrelation.dfa <- lapply(Model_prediction, function(x) mltools::mcc(x, RefT))
    MatthewCorrelation.dfb <- do.call("cbind.data.frame", MatthewCorrelation.dfa)
    MatthewCorrelation.df <- cbind.data.frame(value = as.numeric(MatthewCorrelation.dfb[1,]), L1 = colnames(MatthewCorrelation.dfb) , Statistic.param = "MatthewCorrelation")
    
    ## Adjusted Rank Index for mulitclass
    AdjustedRank.Index.dfa <- lapply(Model_prediction, function(x)  mclust::adjustedRandIndex(x, RefT))
    AdjustedRank.Index.dfb <- do.call("cbind.data.frame", AdjustedRank.Index.dfa)
    AdjustedRank.Index.df <- cbind.data.frame(value = as.numeric(AdjustedRank.Index.dfb[1,]), L1 = colnames(AdjustedRank.Index.dfb), Statistic.param = "AdjustedRank.Index")
    
    ## Table for values
    Global.stats <- rbind.data.frame(Global.stats, MatthewCorrelation.df, AdjustedRank.Index.df)
    colnames(Global.stats) <- c("Value", "Model", "Statistic.parameter")
    KeepStats <- c("Accuracy", "Kappa", "MatthewCorrelation", "AdjustedRank.Index")
    Global.stats <-  Global.stats[Global.stats$Statistic.parameter %in% KeepStats, ] 
    rownames(Global.stats) = NULL
    
    ## List for results
    results = list(Prob.results = res,
                   Confusion.Matrix = Confusion.Matrix,
                   Global.stats = Global.stats,
                   Details.stats = Details.stats,
                   #Discordances.strains = Discordances.strains,
                   #Concordances.strains = Concordances.strains,
                   Correct.ClassificationFreq = Correct.ClassificationFreq,
                   Incorrect.ClassificationFreq = Incorrect.ClassificationFreq)
    
  } else {results = res}
  
  if (is.null(Reference) == TRUE && length(unique(as.character(Reference)))){message("At least two categories are to be compared in peaks object to estimate the statistics")}
  
  return(results)
}