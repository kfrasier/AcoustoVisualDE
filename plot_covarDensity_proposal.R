plot.covarDensity.proposal <- function(mySegments,covarList,presAbs){
  
  # presAbs = a vector the same length as mySegments, 
  #   where 1 = animals were present 0  is no animals were present

  # Figure out which rows had detections, which didn't
  posRows <- which(presAbs==1)
  negRows <- which(presAbs==0)
  
  
  # For each covariate, make density plot with positive and negative encounters
  nCovars <- length(covarList)
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  
  png('density_pres_abs_acoustic_subset.png', width = 2000, height = 800,res=300)
  par(mfrow = c(1,3)) # does the subplots, mar = c(2,3,2,1)
  A <-expression(paste(log[10],'(Mixed Layer Depth) (m)'))
  B <-expression(paste("SST (",~degree*C,")"))
  tList = c(B,A, "Dist. to Eddy (m)");
  i <- 1
  for (cI in covarList){
    if(i==3){
      d0 <- density(mySegments[,cI][negRows]/1000,na.rm =TRUE) # density for negative segments
      d1 <- density(mySegments[,cI][posRows]/1000,na.rm =TRUE) # density for positive segments
    } else {
      d0 <- density(mySegments[,cI][negRows],na.rm =TRUE) # density for negative segments
      d1 <- density(mySegments[,cI][posRows],na.rm =TRUE) # density for positive segments
    }
    x_range = range(c(d0$x,d1$x)) 
    
    if (i==1){
      plot(d1, xlim = x_range, cex.axis = 1,main = "",ylab = "Probabilty of Occurrence", xlab = tList[i])
    }else{
      plot(d1, xlim = x_range, cex.lab = 1,cex.axis = 1,main = "",ylab = "", xlab = tList[i])
    }
    par(new = TRUE)
    plot(d0, xlim = x_range, lty=2,axes = FALSE, xlab = tList[i], ylab = "", main = "")
    
    # Position and make legend
    # legend("topright", c("0","1"), cex = 1,lty = 1:2);
    i <- i+1
  }
  # 
  dev.off()
}