### Predict Category ##

Predict_LogReg <-  function (peaks,
                             model,
                             moz,
                             tolerance = 6,
                             normalizeFun = TRUE,
                             noMatch = 0){

  #Create a subfunction
  predi=function(Peaks,Modell,Moz,Tolerance,NormalizeFun,noMatch){
  ## Empty list
  Peak <- list()
  DF.match <- list()
  DF.Peaks <- list()
  prediction <- list()

  name=rep(NA,length(Peaks))
  method=name

  for (i in 1:length(Peaks)){
    name[i]=Peaks[[i]]@metaData$fullName;
    print(i)
    method[i]= Modell$method;

    Peak[[i]] <- data.frame(X_var = Peaks[[i]]@mass, # mass over charges to match
                            intensity = Peaks[[i]]@intensity) # intensities

    ## Sort the selected mass over charge
    s_moz <- data.frame ("X_var" = unique(sort(Moz)))

    ## Match intensities of the peaks to the selected mass over charge values using the Tolerance argument
    DF.match[[i]] <- fuzzyjoin::distance_left_join(s_moz,
                                                   Peak[[i]],
                                                   by = "X_var",
                                                   max_dist = Tolerance)
    ## If several repetitions of same m/z in selected m/z because of the Tolerance,
    ## only the ones equal in the mass peaks are used if possible
    td=table(DF.match[[i]]$X_var.x)
    dbl=names(td)[td>1]
    to_del=(DF.match[[i]]$X_var.x%in%dbl)*(abs(DF.match[[i]]$X_var.x-DF.match[[i]]$X_var.y)>1e-9)
    to_del[is.na(to_del)]=0
    DF.match[[i]]=DF.match[[i]][(to_del!=1),]

    ## When no match, intensity values are replaced by the noMatch argument
    DF.match[[i]]$intensity[is.na(DF.match[[i]]$intensity)] <- noMatch

    DF.Peaks[[i]] <- data.frame(DF.match[[i]]$intensity)
    DF.Peaks[[i]] <- t(DF.Peaks[[i]])
    colnames(DF.Peaks[[i]]) <- as.character(DF.match[[i]]$X_var.x)

    # Normalize peak with the maximum intensity value
    if (NormalizeFun){
      norma<-function(x) x/(max(x))
      DF.Peaks[[i]] <- apply(DF.Peaks[[i]],1,norma)
      DF.Peaks[[i]] <- t(DF.Peaks[[i]])
    }
    ######
    if (method[i]=="xgbTree"){
      Y_x=DF.Peaks[[i]][,colnames(DF.Peaks[[i]])%in%Modell$finalModel$feature_names]
      Y_new=NULL
      for (k in 1:length(Modell$finalModel$feature_names)){
        Y_new=c(Y_new,Y_x[names(Y_x)==Modell$finalModel$feature_names[k]])
      }
      tY_new=t(Y_new)
      prediction[[i]] <- stats::predict(Modell,tY_new, type = "prob")
    }else{prediction[[i]] <- stats::predict(Modell,DF.Peaks[[i]], type = "prob");}

  }

  ## Results
  Results <- do.call("rbind", prediction)
  Results=cbind(name,method,Results)
  rownames(Results)=paste(1:nrow(Results),".",sep="")

  return(Results)
  }

## Create a data frame for mass spectra Peaks
if (!is.null(model$method)){
   model=list(model)
}
res=NULL
for (i in 1:length(model)){
        res1=predi(Peaks=peaks,Modell=model[[i]],Moz=moz,Tolerance=tolerance,
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
    merge_p=apply(res_1[,3:ncol(res_1)],2,function(x)return(sumlog(as.numeric(x))$p))
    
  #Max vote
    res_1=res[res$name==lr[i],3:ncol(res_1)]
    votes=table(unlist(apply(res_1,1,function(x){return(which(as.numeric(x)==max(as.numeric(x))))})))
    vote_max=which(as.numeric(votes)==max(as.numeric(votes)))
    mat_vote=rep(0,ncol(res_1))
    mat_vote[which(as.numeric(votes)==max(as.numeric(votes)))]=1
    mat_vote=mat_vote/sum(mat_vote)
    
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

return(res)
}
