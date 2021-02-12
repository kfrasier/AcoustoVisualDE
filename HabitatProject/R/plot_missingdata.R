#' Function generates set of histograms showing data gaps by year.
#' RETURNS:
#'   percFilled<- A table of percentage data available by year for each covariate in covarList
#' @export
#' @examples
#' plot.missingdata()

plot.missingdata <- function(mySegments,covarList,fNamePrefix,varUnits){

  # Make vector of dates by year, so we can figure out what % of segments in yeach year are missing data.
  tmpTimes <- as.Date(mySegments$date,"%Y-%m-%d", tz="GMT")
  yearMax <- max(tmpTimes, na.rm=TRUE)
  startYear <- lubridate::floor_date(min(tmpTimes, na.rm=TRUE),unit="year")
  disp(startYear)
  endYear <- lubridate::ceiling_date(max(tmpTimes, na.rm=TRUE),unit="year")
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

  percFilled <- magic::adiag(a,b)


  # for each covariate and year, calculate the percentage of points that have a non-zero value
  m1 <- 1
  # and plot as we go
  nRows = round(sqrt(nCovars))
  nCols = ceiling(nCovars/nRows)
  png(paste0(fNamePrefix,'_missingData.png'), width = 1600, height = 1000)
  par(mfrow=c(nRows,nCols),mai=c(1.2, 1.2, 0.85, 0.42),mgp = c(6, 1, 0))
      # mar = c(5, 6, 4, 2) + 0.1, does the subplots

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
    graphics::barplot(percFilled[m1,]*100, main = varUnits[cI],  xlab = "",
            ylab = "% of data available", las=2, cex.main=3,cex.lab = 3,cex.axis = 2.5,cex.names = 2.5,
            ylim=c(0,100))
    m1 <- m1 + 1
  }

  dev.off()
  # return(percFilled)

}
