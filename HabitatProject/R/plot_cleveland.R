#' Function plots covariates to look for skew 
#' @export 
#' @examples 
#' plot.cleveland() 

plot.cleveland <- function(mySegments,covarList,isTrans,fNamePrefix,varUnits){
  # isTrans is a boolean flag to name the file appropriately
  
  nCovars <- length(covarList)
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  

  if (!isTrans){
    png(paste0(fNamePrefix,"_clevelandDots_noTransform.png"), width = 1000, height =1000)
  }else{
    png(paste0(fNamePrefix,"_clevelandDots_withTransform.png"), width = 1000, height =1000)
    
  }
  
  par(mfrow = c(nRows,nCols))#,mai=c(.6, .6, 0.6, 0.42))# mar = c(2,1,2,1)) # does the subplots
  
  for (cI in covarList){
    thisCovar <- mySegments[,cI]
    dotchart(thisCovar, xlab = varUnits[cI], cex.lab = 1.5, cex = 1)
  }
  dev.off()
  
  
}
