
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

    if (length(unique(names(peaks))) != length(peaks) && is.null(names(peaks)) == FALSE) {
      stop("Each element of X must have a unique name !")
    }

    predi = function(Peaks, Modell, Moz, Tolerance, NormalizeFun, noMatch) {
      Peak <- list()
      DF.match <- list()
      DF.Peaks <- list()
      prediction <- list()
      diff_mass <- list()
      DF.matchSplit <- list()
      name = rep(NA, length(Peaks))
      method = name
      Tolerance.step <- list()
      for (i in 1:length(Peaks)) {
        if (is.null(Peaks[[i]]@metaData$fullName)) {
          name[i] = Peaks[[i]]@metaData$file
        }
        else {
          name[i] = Peaks[[i]]@metaData$fullName
        }
        print(name[i])
        method[i] = Modell$method
        Peak[[i]] <- data.frame(X_var = Peaks[[i]]@mass,intensity = Peaks[[i]]@intensity)
        s_moz <- data.frame(X_var = unique(sort(as.numeric(Moz))))
        DF.match[[i]] <- d_left_join(s_moz,Peak[[i]], by = "X_var", max_dist = Tolerance)
        if (sum(is.na(DF.match[[i]]$X_var.y)) == nrow(DF.match[[i]])) {
          warning("No m/z from peaks object matched with the m/z in the moz objet,\n
                  the tolerance has been increased by steps indicated in toleranceStep argument Da")
          Tolerance.step <- Tolerance
          while (sum(is.na(DF.match[[i]]$X_var.y)) == nrow(DF.match[[i]])) {
            Tolerance.step <- Tolerance.step + toleranceStep
            DF.match[[i]] <- d_left_join(s_moz,Peak[[i]], by = "X_var", max_dist = Tolerance.step)
            Tolerance.step <- Tolerance.step
          }
          message(paste(c("tolerance found for match =",Tolerance.step)))
        }
        DF.match[[i]]$intensity[is.na(DF.match[[i]]$intensity)] <- noMatch
        td = table(DF.match[[i]]$X_var.x)
        dbl = names(td)[td > 1]
        if (length(dbl) > 0) {
          keepDF = NULL
          for (l in 1:length(dbl)) {
            diff_m = abs(DF.match[[i]][DF.match[[i]]$X_var.x == dbl[l], 2] - DF.match[[i]][DF.match[[i]]$X_var.x == dbl[l], 1])
            keepDF = rbind(keepDF, DF.match[[i]][DF.match[[i]]$X_var.x == dbl[l], ][which.min(diff_m), ])
          }
          DF.match[[i]] = DF.match[[i]][which(!DF.match[[i]]$X_var.x %in% dbl), ]
          DF.match[[i]] = rbind(DF.match[[i]], keepDF)
          DF.match[[i]] = DF.match[[i]][order(DF.match[[i]][,1]), ]
        }
        DF.Peaks[[i]] <- data.frame(DF.match[[i]]$intensity)
        DF.Peaks[[i]] <- t(DF.Peaks[[i]])
        colnames(DF.Peaks[[i]]) <- as.character(DF.match[[i]]$X_var.x)
        if (NormalizeFun) {
          norma <- function(x) x/(max(x))
          DF.Peaks[[i]] <- apply(DF.Peaks[[i]], 1, norma)
          DF.Peaks[[i]] <- t(DF.Peaks[[i]])
        }
        prediction[[i]] <- stats::predict(Modell, DF.Peaks[[i]], type = "prob")
      }
      Results <- do.call("rbind", prediction)
      Results = cbind(name, method, Results)
      rownames(Results) = paste(1:nrow(Results), ".",
                                sep = "")
      NameStrain <- rep(names(peaks), nlevels(factor(method)))
      if (is.null(names(peaks))) {
        Results$name = rep(c(1:length(peaks)), nlevels(factor(method)))
      }
      else {
        Results$name = rep(NameStrain, nlevels(factor(method)))
      }
      return(Results)
    }
    if (!is.null(model$method)) {
      model = list(model)
    }
    res = NULL
    for (i in 1:length(model)) {
      res1 = predi(Peaks = peaks, Modell = model[[i]], Moz = as.numeric(moz),
                   Tolerance = tolerance, NormalizeFun = normalizeFun,
                   noMatch = noMatch)
      res1$method=paste(res1$method,i,sep="_")
      res = rbind(res, res1)
    }
    if (length(model) > 1) {
      lr = levels(as.factor(res$name))
      if (length(lr) > 1) {
        for (i in 1:length(lr)) {
          res_1 = res[res$name == lr[i], ]
          suppressWarnings({
            merge_p = apply(res_1[, 3:ncol(res_1)], 2, function(x) return(metap::sumlog(as.numeric(x))$p))
          })
          res_1 = res[res$name == lr[i], 3:ncol(res_1)]
          Max1 = apply(data.frame(res_1), 1, function(x) replace(as.numeric(x),which.max(x), 1))
          Max2 = apply(Max1, 2, function(x) replace(x,x < 1, 0))
          VoteTot = ncol(Max2)
          Votedecompte = apply(Max2, 1, sum)
          Votedecompte2 = Votedecompte/VoteTot
          mat_vote = Votedecompte2
          res = rbind(res, c(lr[i], "comb_fisher",merge_p))
          res = rbind(res, c(lr[i], "max_vote", mat_vote))
        }
      }
    }
    pred_cat = NULL
    for (i in 1:nrow(res)) {
      pred_cat = c(pred_cat, colnames(res)[3:ncol(res)][which.max(res[i,3:ncol(res)])])
    }
    res = cbind(res, pred_cat)
    colnames(res)[ncol(res)] = "pred_max_p"
    if (is.null(Reference) == FALSE && length(unique(as.character(Reference))) >  1) {
      df_methods <- split.data.frame(res, res$method)
      suppressWarnings({
        RefT = factor(make.names(Reference))
        Model_prediction <- lapply(df_methods, function(x) factor((x[, "pred_max_p"]), levels = levels(RefT)))
        Confusion.Matrix <- lapply(Model_prediction, function(x) try(caret::confusionMatrix(x, RefT), silent = TRUE))
        Confusion.MatrixT <- lapply(Confusion.Matrix, function(x) if (is.character(x) == FALSE) {
          x <- x
        })
        Confusion.MatrixT <- Confusion.MatrixT[!sapply(Confusion.MatrixT,is.null)]
        Confusion.Matrix = Confusion.MatrixT
      })
      Df.classif <- lapply(df_methods, function(x) try(cbind.data.frame(Reference = make.names(Reference), x), silent = TRUE))
      Df.classif <- Df.classif[!sapply(Df.classif, is.character)]
      DiscordTable <- lapply(Df.classif, function(x) x[which(x[, "Reference"] != x[, "pred_max_p"]), ])
      ConfMat.df <- lapply(Confusion.Matrix, function(x) data.frame(x[["table"]]))
      Discord2 <- lapply(ConfMat.df, function(x) x[which(x[,1] != x[, 2]), ])
      Incor.Classification <- lapply(Discord2, function(x) x[which(x[,3] >= 1), ])
      suppressMessages(Incorrect.ClassificationFreq <- try(reshape2::melt(Incor.Classification)[,-3], silent = TRUE))
      if (is.character(Incorrect.ClassificationFreq)) {
        Incorrect.ClassificationFreq <- Incorrect.ClassificationFreq[!sapply(Incorrect.ClassificationFreq,is.character)]
      }
      if (is.character(Incorrect.ClassificationFreq)) {
        Incorrect.ClassificationFreq <- t(data.frame(c(0,0, 0, 0)))
      }
      colnames(Incorrect.ClassificationFreq) <- c("Prediction","Reference", "Frequence", "Model")
      ConcordTable <- lapply(Df.classif, function(x) x[which(x[,"Reference"] == x[, "pred_max_p"]), ])
      Cor.Classification <- lapply(ConfMat.df, function(x) x[which(x[,1] == x[, 2]), ])
      Cor.Classification <- lapply(Cor.Classification, function(x) x[which(x[,3] > 0), ])
      suppressMessages(Correct.ClassificationFreq <- reshape2::melt(Cor.Classification)[,-3])
      colnames(Correct.ClassificationFreq) <- c("Prediction", "Reference", "Frequence", "Model")
      Global.stat <- lapply(Confusion.Matrix, function(x) (x[["overall"]]))
      Global.stats <- do.call("data.frame", Global.stat)
      Statistic.param = row.names(Global.stats)
      Global.stats <- suppressMessages(reshape2::melt(Global.stat,id.vars = NULL))
      Global.stats <- data.frame(Global.stats, Statistic.param)
      Details.stat <- try(lapply(Confusion.Matrix, function(x) cbind.data.frame(x[["byClass"]],Class = row.names(x[["byClass"]]))), silent = TRUE)
      suppressMessages(Details.stats <- try(reshape2::melt(Details.stat),silent = TRUE))
      if (is.character(Details.stat)) {
        Details.stat2 <- lapply(Confusion.Matrix, function(x) cbind.data.frame(x[["byClass"]],Class = x[["positive"]]))
        suppressMessages(Details.stats2 <- cbind.data.frame(reshape2::melt(Details.stat2)))
        colnames(Details.stats2) <- c("Class", "Statistic.parameter","Value", "Model")
        Details.stats2$Statistic.parameter <- row.names(data.frame(Details.stat2[1]))
        Details.stats = Details.stats2
      }
      MatthewCorrelation.dfa <- lapply(Model_prediction, function(x) mltools::mcc(x,RefT))
      MatthewCorrelation.dfb <- do.call("cbind.data.frame",MatthewCorrelation.dfa)
      MatthewCorrelation.df <- cbind.data.frame(value = as.numeric(MatthewCorrelation.dfb[1,]),
                                                L1 = colnames(MatthewCorrelation.dfb), Statistic.param = "MatthewCorrelation")
      AdjustedRank.Index.dfa <- lapply(Model_prediction, function(x){mclust::adjustedRandIndex(x,RefT);})
      AdjustedRank.Index.dfb <- do.call("cbind.data.frame",AdjustedRank.Index.dfa)
      AdjustedRank.Index.df <- cbind.data.frame(value = as.numeric(AdjustedRank.Index.dfb[1,]),
                                                L1 = colnames(AdjustedRank.Index.dfb), Statistic.param = "AdjustedRank.Index")
      Global.stats <- rbind.data.frame(Global.stats, MatthewCorrelation.df,AdjustedRank.Index.df)
      colnames(Global.stats) <- c("Value", "Model","Statistic.parameter")
      KeepStats <- c("Accuracy", "Kappa", "MatthewCorrelation","AdjustedRank.Index")
      Global.stats <- Global.stats[Global.stats$Statistic.parameter %in% KeepStats, ]
      rownames(Global.stats) = NULL
      results = list(Prob.results = res, Confusion.Matrix = Confusion.Matrix,
                     Global.stats = Global.stats, Details.stats = Details.stats,
                     Correct.ClassificationFreq = Correct.ClassificationFreq,
                     Incorrect.ClassificationFreq = Incorrect.ClassificationFreq)
    }
    else {
      results = res
    }
    if (is.null(Reference) == TRUE && length(unique(as.character(Reference)))) {
      message("At least two categories have to be compared in peaks objects to estimate the statistics.")
    }
  return(results)
}
