#' Function plots presence absence as a function of predictor variables. 
#' @export 
#' @examples 
#' plot.covarDensity() 

plot.covarDensity <- function(mySegments,covarList,presAbs,fNamePrefix,varUnits){
  
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
  par(mfrow = c(nRows,nCols), mar = c(4,4,2,1)) # does the subplots
  
  
  i <- 1
  for (cI in covarList){
    d0 <- stats::density(mySegments[,cI][negRows],na.rm = TRUE) # density for negative segments
    d1 <- stats::density(mySegments[,cI][posRows],na.rm = TRUE) # density for positive segments
    
    x_range = range(c(d0$x,d1$x)) 
    
    graphics::plot(d0, xlim = x_range, cex.main = 1.5, cex.lab=1.5,
                     cex.axis = 1.5, main = "", xlab = varUnits[cI], ylab = 'normalized counts')
    par(new = TRUE)
    graphics::plot(d1, xlim = x_range, lty=2,axes = FALSE, xlab = "", ylab = "", main = "")
    
    # Position and make legend
    legend("topright", c("0","1"), cex = 1.5,lty = 1:2,bty = "n");
    i <- i+1
  }
  # 
  dev.off()
}