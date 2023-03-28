########################
### SelectionVarStat ##
########################

SelectionVarStat=function(X,
                          Y,
                          stat.test="Limma",
                          pi0.method="abh",
                          fdr=0.05,
                          Sampling = c("no", "up","down", "smote")){
  
  
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
         })
  
  
  
  #Non parametric test 
  if (stat.test=="kruskal"){
    vp=apply(X,2,function(x)stats::kruskal.test(x,g = Y)$p.value)
  }
  
  #Contrast matrix for parametric tests
  ct=matrix(0,length(levels(Y)),length(levels(Y))-1)
  ct[2,1]=1
  if (nrow(ct)>=3){
    for (j in 2:(nrow(ct)-1)){
      ct[j,j]=1
      ct[j+1,j]=-1
    }}
  
  #ANOVA test
  if (stat.test=="anova"){
    p.value_ANOVA=rep(NA, ncol(X))
    for (i in 1:ncol(X)){
      mod <- stats::lm(X[,i] ~ Y)
      p.value_ANOVA[i]=car::linearHypothesis(mod,t(ct))$'Pr(>F)'[2]
    } 
    vp=p.value_ANOVA
  }
  
  #Limma test
  if (stat.test=="Limma"){
    design=stats::model.matrix(~Y)
    fit=limma::eBayes(limma::contrasts.fit(limma::lmFit(t(X), design),ct),robust=TRUE)
    p.value_LIMMA=fit$F.p.value 
    vp=p.value_LIMMA
  }
  
  #Estimate proportion of "not moving" features
  pi0=cp4p::estim.pi0(vp,pi0.method = pi0.method)
  
  #Expected nb of features to select
  nb_to_sel=floor(ncol(X)*(1-pi0[[1]]))
  
  #P-value adjustment
  ap=cp4p::adjust.p(vp,pi0.method=pi0.method)
  
  #Selected features using the fdr threshold
  sel_moz=colnames(X)[which(ap$adjp$adjusted.p<fdr)]
  
  return(list(nb_to_sel=nb_to_sel,sel_moz=sel_moz))
}