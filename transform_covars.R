# plot covariates to look for skew
transform.covars <- function(mySegments,covarList,transformList){
  # Transform list must be the same length as the covariate list!
  dataNames <- NULL
  
  tI <- 1
  for (cI in covarList){
    
    thisCovar <- mySegments[,cI]
    thisTransform <- transformList[tI]
    
    if (thisTransform == "log10"){
      thisCovarTr <- log10(thisCovar)
      dataNames[tI] <- paste0('log10_',cI,sep = '')
    }else if (thisTransform == "ln"){
      thisCovarTr <- log(thisCovar)
      dataNames[tI] <- paste0('ln_',cI,sep = '')
    }else if (thisTransform == "exp"){
      thisCovarTr <- exp(thisCovar)
      dataNames[tI] <- paste0('exp_',cI,sep = '')
    }else if (thisTransform == "sqrt"){
      thisCovarTr <- sqrt(thisCovar)
      dataNames[tI] <- paste0('sqrt_',cI,sep = '')
    }else if (thisTransform == "^2"){
      thisCovarTr <- (thisCovar^2)
      dataNames[tI] <- paste0('sq_',cI,sep = '')
    }else if (thisTransform == "none"){
      thisCovarTr <- thisCovar
      dataNames[tI] <- cI
    }else{
      cat(paste("Transform on ", cI, " not recognized.\n"))
      thisCovarTr <- thisCovar
      dataNames[tI] <- cI
    }
    
    if (tI == 1){
      covarTransform <- data.frame(thisCovarTr)
    }  else{ 
      covarTransform[,cI] <- thisCovarTr  
    }
 
    tI <- tI + 1
  }
  colnames(covarTransform) <- dataNames
  
  return(covarTransform)
  
}