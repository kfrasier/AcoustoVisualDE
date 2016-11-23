# fit_visdata
fit.visdata<- function(visData,spIdxON2,PLC,keyList1,adjList1){

  # Processes Visual data to determine detection probabilities and strip widths for each platform
  # Saves model output to files, and makes plots for evaluating model fit.
  # Returns: 
  # tDist <- Truncation distances
   
  # for each visual platform
  nPlatform <- 1  # Initialize to start with first platform
  
  tDist <-NULL # store truncation distances
  bestModel <-NULL # store best model covariates
  bestKey <-NULL # store best model key
  
  cat("Begin model fitting for visual data \n")
  for (i in PLC){
    ddfData <- NULL
    noquote(paste("Fitting platform ", i))
    
    # identify sightings assoicated with the a certain platform
    PLCspIdxOn <- spIdxON2[which(grepl(i,visData$ship[spIdxON2]))] 
    
    # Make dataframe with the inputs that the ddf distance function wants
    ddfData$observer <- rep(1,length(PLCspIdxOn))
    ddfData$detected <- rep(1,length(PLCspIdxOn))
    ddfData$object <- (1:length(PLCspIdxOn))
    ddfData$distance <- visData$transect_distm[PLCspIdxOn]
    ddfData$size <- visData$size[PLCspIdxOn]
    ddfData$seastate <- visData$seastate[PLCspIdxOn]
    ddfData$swell <- visData$swell[PLCspIdxOn]
    ddfData$vis <- visData$vis[PLCspIdxOn]
    
    
    # Compute untruncated detection function
    cat("Calculating basic fit with non-truncated data and half-normal key, no covariates.\n")
    
    detFun_noTrunc <-ddf(method='ds',dsmodel=~mcds(key='hn', formula = ~ 1),
                         data=as.data.frame(ddfData), meta.data=list(binned=F, left=0))
    
    
    # Make output plots
    noquote("Saving plots")
    png(paste(SP,'sightnoTrunc_',i,'.png',sep=''), width = 800, height = 500)
    par(mfrow=c(1,2))
    plot(detFun_noTrunc)
    qqplot.ddf(detFun_noTrunc,plot=TRUE)
    dev.off()
    
    
    # Compute truncation distance by removing highest 5% of distances
    tDist[nPlatform] <- quantile(ddfData$distance,.95,na.rm = TRUE)
    cat(paste("Truncation distance for platform ", i, "=",  round(tDist[nPlatform],2), "km \n"))
    
    
    # Iterate over detection functions with various adjustments and orders, and identify AIC for each
    # list of key funs to try: 
    # keyList1 = list(  "hn", "hn",  "hr",  "hr",  "hr","unif","unif")
    # adj  = c("none","cos","none","poly","poly","none","cos")
    adjOrder = c(     0,    3,     0,     2,     4,     0,    2)
    adjStr = paste(adjList1,adjOrder)
    detFun1 <- NULL
    aicList1 <- NULL
    
    cat("Fitting detection functions with adjustments \n")
    
    for (i1 in 1:length(keyList1)){
      if (grepl('none',adjList1[[i1]])){
        
      
        detFun1[[i1]] <-ddf(method ='ds', dsmodel =~ mcds(key = as.character(keyList1[i1]), formula = ~ 1),
                            data = as.data.frame(ddfData), meta.data = list(binned=F, width=tDist[nPlatform], left=0))
        
      } else {
        detFun1[[i1]]<-ddf(method ='ds', dsmodel =~ mcds(key = as.character(keyList1[i1]), formula = ~ 1,
                           adj.series = as.character(adjList1[i1]), adj.order = adjOrder[i1]), data = as.data.frame(ddfData),
                           meta.data = list(binned=F, width=tDist[nPlatform],left=0))
      }
      aicList1[i1] <- detFun1[[i1]]$criterion
      cat(paste0("Model result ", i1,": Key = ",keyList1[i1],
                 "; Adjustment = ",adj[i1], "; Order = ",adjOrder[i1],"\n", collapse = ""))
      cat(paste("AIC =",  round(aicList1[i1], digits=2),"\n"))
      
    }
    # Iterate over detection functions with covariates and identify AIC for each
    covarList <- c('size', 'seastate', 'vis', 'swell')
    
    # list of keys: hn, hr
    keyCovar= c('hn', 'hr')
    CI <- 1
    aicList2 <- NULL # store AIC scores 
    keyList2 <- NULL # store keys scores 
    cSetStr <- NULL # store the covariate formulas
    detFun2<- NULL
    # iterate over the key options
    for (iKey in keyCovar){
      
      # iterate over covariate combinations - using "combn" to come up with different the combinations of covariates
      for (i2 in 1 : length(covarList)){
        covarSet <- combn(covarList,i2)
        nSets <- dim(covarSet)
        
        for (i3 in 1:nSets[2]){
          cSet <- covarSet[,i3]
          # nPlus <- nSets[1]-1
          cSetStr[CI] <- paste('~', paste0(cSet, collapse = " + "))
          
          # sometimes models do not converge, use try() to avoid crash if a model fails
          dF <- NULL
          try(dF <- ddf(method='ds',dsmodel=~mcds(key=iKey,  formula = cSetStr[CI]),
                        data = as.data.frame(ddfData),
                        meta.data = list(binned=F, width=tDist[nPlatform],left=0)))
          
          if (is.null(dF)){
            cat(paste0("Model did not converge: Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
            
          }else {
            detFun2[[CI]] <-dF
            
            aicList2[CI] <- detFun2[[CI]]$criterion
            keyList2[CI] <- iKey
            cat(paste0("Model result ", CI+i1, ": Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
            cat(paste("AIC =",  round(aicList2[CI], digits=2)))
          }
          
          CI = CI+1
          
        }
      }
    }
    cat("Done fitting models")
    
    # Put all combinations together, and see which one has the lowest AIC
    aicList<-c(aicList1,aicList2)
    keyList <- c(keyList1,keyList2)
    adjList <- c(adjStr,cSetStr)
    detFun <- c(detFun1,detFun2)
    
    ddfOut <- data.frame(model = adjList, key = keyList, aic = aicList)
    
    # Best model is...
    bestModelIdx <- which(ddfOut$aic == min(ddfOut$aic, na.rm = TRUE), arr.ind = TRUE)
    bestModel[[nPlatform]] <- adjList[bestModelIdx]
    bestKey[[nPlatform]] <- keyList[bestModelIdx]
    cat(paste("Best model for Platform ", i,":\n"))
    cat(paste("Key = ", bestKey[[nPlatform]], "; Adjustment =",  bestModel[[nPlatform]],"\n"))
    
    
    # Make output plot of best model
    cat("Saving plots and summaries \n")
    png(paste(SP,'sightwTrunc_',i,'.png',sep=''), width = 800, height = 500)
    par(mfrow=c(1,2))
    plot(detFun[[bestModelIdx]], main = paste('model = ',bestModel[[nPlatform]], '; key = ', bestKey[[nPlatform]]))
    qqplot.ddf(detFun[[bestModelIdx]],plot=TRUE)
    dev.off()
    
    # Output summary text to txt file
    sink(paste(SP,'sightwTrunc_',i,'.txt',sep=''))
    summary(detFun[[bestModelIdx]])
    ddf.gof(detFun[[bestModelIdx]])
    sink()
    
    ## save best model
    save(detFun, bestModelIdx, file = paste(SP,'sightwTrunc_',i,'.Rdata',sep=''))
    
    nPlatform <- nPlatform + 1
    cat(paste("Done fitting models for platform ", i,"\n"))
  }
  
  cat("Done fitting models \n")
  
  
  return(tDist)
}


