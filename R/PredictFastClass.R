

PredictFastClass <- function(peaks, 
                             mod_peaks,
                             Y_mod_peaks,
                             moz="ALL",
                             tolerance = 6,
                             toleranceStep = 2,
                             normalizeFun = TRUE,
                             noMatch = 0){
  
  if (isTRUE(moz=="ALL")){moz=colnames(mod_peaks)}else{moz=moz}
  if (length(levels(as.factor(as.character(Y_mod_peaks))))==1){warning("There is only one category in the set of mass spectra used in regression models !");}
  if (nrow(mod_peaks)!=length(Y_mod_peaks)){warning("Number of MS in mod_peaks does not correspond to the length of Y_mod_peaks !");}
  
  Peaks=peaks
  ## Empty list
  Peak <- list()
  DF.match <- list()
  DF.Peaks <- list()
  prediction <- list()
  diff_mass <- list()
  DF.matchSplit <- list()
  Tolerance.step <- list()
  
  
  p.value_cat=NULL
  p.value_cat2=NULL

  nam=rep(0,length(Peaks))
  
  for (i in 1:length(Peaks)){

    Peak[[i]] <- data.frame(X_var = Peaks[[i]]@mass, # mass over charges to match
                            intensity = Peaks[[i]]@intensity) # intensities
    
    ## Sort the selected mass over charge
    s_moz <- data.frame ("X_var" = unique(sort(as.numeric(moz))))
    
    ## Match intensities of the peaks to the selected mass over charge values using the Tolerance argument
    DF.match[[i]] <- fuzzyjoin::distance_left_join(s_moz,
                                                   Peak[[i]],
                                                   by = "X_var",
                                                   max_dist = tolerance)
    
    ### If no match
    if(isTRUE(sum(DF.match[[i]]$X_var.y, na.rm = TRUE) == 0)){
      warning("No m/z from peaks object matched with the m/z in the moz objet,
                the tolerance has been increased by steps indicated in toleranceStep argument Da")
      
      Tolerance.step[[i]] <- tolerance + toleranceStep  #Tolerance step
      
      while(isTRUE(sum(DF.match[[i]]$X_var.y, na.rm = TRUE) == 0)){
        
        
        DF.match[[i]] <- fuzzyjoin::distance_left_join(s_moz,
                                                       Peak[[i]],
                                                       by = "X_var",
                                                       max_dist = Tolerance.step[[i]])
        
        Tolerance.step[[i]] <- Tolerance.step[[i]] + toleranceStep #Tolerance step
      }
      message(paste(c("tolerance found for match=", Tolerance.step[[i]])))
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
      if (normalizeFun){
        norma<-function(x) x/(max(x))
        DF.Peaks[[i]] <- apply(DF.Peaks[[i]],1,norma)
        DF.Peaks[[i]] <- t(DF.Peaks[[i]])
      }

    
  
    new_peaks=t(DF.Peaks[[i]])
    nam[i]=Peaks[[i]]@metaData$fullName
    # Print Metadata
    print(nam[i])
    IntM=t(mod_peaks[,which(colnames(mod_peaks)%in%s_moz$X_var)])
    Y=Y_mod_peaks
  
    
    newXIC=as.matrix(data.frame(new_peaks))[,1]
    p.value_cat=rbind(p.value_cat,rep(NA,length(levels(Y))))
    p.value_cat2=rbind(p.value_cat2,rep(NA,length(levels(Y))))
    if (length(newXIC)!=nrow(IntM)){
      warning("Problem: Mass Spectrum and MS dataset does not have same length of m/z.")
    }else{
      sna=apply(IntM[which(!is.na(newXIC)),],1,function(x){sum(is.na(x));})
      if (length(which(sna!=ncol(IntM)))==0){
        warning("Mass Spectrum does not match any mass spectra in the MS dataset: it is probably from another category.")
        p.value_cat2[k,]=rep(1,length(levels(Y)))
        p.value_cat[k,]=rep(NA,length(levels(Y)))
      }else{
        #test si on a une correspondance parfaite
        se=apply(IntM,2,function(x){sum((x-newXIC)^2,na.rm=T);})
        nbna=apply(IntM,2,function(x){sum(is.na(x-newXIC));})
        se2=se[nbna!=length(newXIC)];
        if (length(which(se2==0))>0){
          warning(paste(c("Mass Spectrum is matching perfectly MS in the database")));
          p.value_cat2[i,]=rep(1,length(levels(Y)))
          p.value_cat2[i,which(levels(Y)==Y[which((se==0)&(nbna!=length(newXIC)))])]=0;
          p.value_cat[i,]=rep(0,length(levels(Y)))
          p.value_cat[i,which(levels(Y)==Y[which((se==0)&(nbna!=length(newXIC)))])]=-1;
        }else{
          desig = stats::model.matrix(~Y-1)
          for (k in 1:ncol(desig)){
            if (length(which(desig[,k]==1))>1){
            sna=apply(IntM[which(!is.na(newXIC)),which(desig[,k]==1)],1,function(x){sum(is.na(x));})
            }else{
              sna=sum(is.na(IntM[which(!is.na(newXIC)),which(desig[,k]==1)]));
              warning(paste(c("You do not have enough spectra in a category ! (at least 2 are required)")));
            }
            if (length(which(sna!=ncol(IntM[,which(desig[,k]==1)])))==0){
              p.value_cat2[i,k]=1;
              p.value_cat[i,k]=NA
            }else{
              newXIC2=newXIC[!is.na(newXIC)]
              B=IntM[!is.na(newXIC),which(desig[,k]==1)]
              B[is.na(B)]=0
              mod = lm(I(newXIC2*1e10) ~ B)
              if (mod$df.residual==0){
                p.value_cat2[i,k]=NA
                p.value_cat[i,k]=NA
              }else{
                if (sum(is.na(mod$coefficients))>0){
                  mod_wo_aliased=stats::lm(I(newXIC2*1e10) ~ B[,which(!is.na(mod$coefficients))-1])
                }else{
                  mod_wo_aliased=mod;
                }
                #F test
                smod=summary(mod_wo_aliased)
                p.value_cat2[i,k]= stats::pf(smod$fstatistic[1],smod$fstatistic[2],smod$fstatistic[3],lower.tail=FALSE)
                #AIC
                p.value_cat[i,k]= stats::AIC(mod_wo_aliased)
                
              }
            }
          }
        }
      }
    }
  }
  
  colnames(p.value_cat)=levels(Y)
  colnames(p.value_cat2)=levels(Y)
  res=data.frame(nam,p.value_cat)
  #Add a column indicating the predicted category and the minimum p-value
  pred_cat=NULL
  min_p=NULL
  for (i in 1:nrow(res)){
    if (length(which.min(res[i,2:ncol(res)]))!=0){
      pred_cat=c(pred_cat,colnames(res)[2:ncol(res)][which.min(res[i,2:ncol(res)])]);
    }else{
      pred_cat=c(pred_cat,NA);
    }
    min_p=c(min_p,min(p.value_cat2[i,]))
  }
  res=cbind(res,min_p,pred_cat)
  
  colnames(res)[1]="name"
  colnames(res)[ncol(res)-1]="p_not_in_DB"
  colnames(res)[ncol(res)]="pred_cat"
  return(res)
}



