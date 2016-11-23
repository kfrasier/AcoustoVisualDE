# plot covariates to look for skew
plot.cleveland <- function(mySegments,covarList,isTrans){
  # isTrans is a boolean flag to name the file appropriately
  
  nCovars <- length(covarList)
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  

  if (!isTrans){
    png("clevelandDots_noTransform.png")
  }else{
    png("clevelandDots_withTransform.png")
    
  }
  
  par(mfrow = c(nRows,nCols),mar = c(2,1,2,1)) # does the subplots
  
  for (cI in covarList){
    thisCovar <- mySegments[,cI]
    dotchart(thisCovar, main = cI)
  }
  dev.off()
  
  
}
