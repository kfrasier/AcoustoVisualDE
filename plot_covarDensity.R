plot.covarDensity <- function(mySegments,covarList,presAbs,fNamePrefix){
  
  # presAbs = a vector the same length as mySegments, 
  #   where 1 = animals were present 0  is no animals were present

  # Figure out which rows had detections, which didn't
  posRows <- which(presAbs>0)
  negRows <- which(presAbs==0)
  
  
  # For each covariate, make density plot with positive and negative encounters
  nCovars <- length(covarList)
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  
  png(paste0(fNamePrefix,'_density_pres_abs.png'), width = 1000, height = 600)
  par(mfrow = c(nRows,nCols), mar = c(2,3,2,1)) # does the subplots
  
  
  i <- 1
  for (cI in covarList){
    d0 <- density(mySegments[,cI][negRows],na.rm = TRUE) # density for negative segments
    d1 <- density(mySegments[,cI][posRows],na.rm = TRUE) # density for positive segments
    
    x_range = range(c(d0$x,d1$x)) 
    
    plot(d0, xlim = x_range, main = cI, cex.main = 1.5, cex.axis = 1.5)
    par(new = TRUE)
    plot(d1, xlim = x_range, lty=2,axes = FALSE, xlab = "", ylab = "", main = "")
    
    # Position and make legend
    legend("topright", c("0","1"), cex = 1.5,lty = 1:2,bty = "n");
    i <- i+1
  }
  # 
  dev.off()
}