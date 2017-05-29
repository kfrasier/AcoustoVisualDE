plot.missingdata <- function(mySegments,covarList,fNamePrefix){
  # Generates set of histograms showing data gaps by year
  # RETURNS:
  #  percFilled<- A table of % data available by year for each covariate in covarList
  
  # Make vector of dates by year, so we can figure out what % of segments in yeach year are missing data.
  tmpTimes <- as.Date(mySegments$date,"%Y-%m-%d")
  yearMax <- max(tmpTimes)
  startYear <- floor_date(min(tmpTimes),unit="year")
  endYear <- ceiling_date(max(tmpTimes),unit="year")
  yearVector <- seq(startYear, endYear, by = "years")
  
  
  # Figure out which year goes with each row of Segments
  yearIdx <- findInterval(tmpTimes,yearVector)
  uniqueYears <- sort(unique(yearIdx))
  
  nYears <- length(uniqueYears)
  nCovars <- length(covarList)
  # Make a nice matrix to store these percentages in
  a <- array(0,dim=c(0,nYears))
  dimnames(a) <- list(covariate = c(),year = as.character(format(yearVector[uniqueYears],"%Y")))
  b <- array(0,dim=c(nCovars,0))
  dimnames(b) <- list(covariate = covarList,year = c())
  
  percFilled <- adiag(a,b)
  
  
  # for each covariate and year, calculate the percentage of points that have a non-zero value
  m1 <- 1
  # and plot as we go
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  png(paste0(fNamePrefix,'_missingData.png'), width = 1600, height = 1000)
  par(mfrow=c(nRows,nCols)) # does the subplots
  
  for (cI in covarList){
    thisCovar <- mySegments[,cI]
    
    n1 <- 1
    for (yI in uniqueYears){
      # Split by year, and ask what % of points in this year have a value
      thisYearSet <- which(yearIdx == yI)
      thisYearVals <- thisCovar[thisYearSet]
      # nFilled <-which(thisYearVals!=-99999)
      nFilled <- which(!is.na(thisYearVals))
      percFilled[m1,n1] <- length(nFilled)/length(thisYearVals)
      n1 <- n1 + 1
      
    }
    barplot(percFilled[m1,]*100, main = cI,  xlab = "",
            ylab = "% of data available", las=2, cex.main=1.5,cex.axis = 1.5,cex.names = 1.5,
            ylim=c(0,100))
    
    m1 <- m1 + 1
  }
  
  dev.off()
  return(percFilled)
  
}