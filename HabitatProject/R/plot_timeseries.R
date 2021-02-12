#' Function to plot acoustic timeseries data and save to image.
#' @export
#' @examples
#' plot.timeseries()

plot.timeseries <-function(siteList,outDir,AcOnlySegments){
  deplGapsStart<-NULL
  deplGapsEnd<-NULL

  deplGapsStart$MC <- as.POSIXct(as.Date(c( "2011-8-14","2013-8-04")),tz = "GMT")

  deplGapsEnd$MC <- as.POSIXct(as.Date(c("2011-9-22","2014-01-1")),tz = "GMT")

  deplGapsStart$GC <- as.POSIXct(as.Date(c("2011-2-3","2011-8-8","2013-9-11")),tz = "GMT")

  deplGapsEnd$GC <- as.POSIXct(as.Date(c("2011-3-22","2011-9-23","2014-1-1")),tz = "GMT")

  deplGapsStart$DT <- as.POSIXct(as.Date(c("2011-1-1","2011-6-25","2011-11-15","2012-1-10","2013-8-19")),tz = "GMT")

  deplGapsEnd$DT <- as.POSIXct(as.Date(c("2011-3-3","2011-7-10","2011-12-13","2012-5-27","2013-10-31")),tz = "GMT")

  deplGapsStart$DC <- as.POSIXct(as.Date(c("2011-02-07","2011-7-7","2012-3-3","2012-12-10","2013-9-26")),tz = "GMT")

  deplGapsEnd$DC <- as.POSIXct(as.Date(c("2011-3-21","2011-10-25","2012-3-3","2012-12-9","2013-12-18")),tz = "GMT")

  deplGapsStart$MP <- as.POSIXct(as.Date(c("2011-02-20","2011-09-07","2012-03-02","2012-10-21","2013-09-26")),tz = "GMT")

  deplGapsEnd$MP <- as.POSIXct(as.Date(c("2011-05-20","2011-09-22","2012-02-29","2012 12 10","2014-1-1")),tz = "GMT")


  for (iSite in 1:length(siteList)){
    if (iSite == length(siteList)){
      png(paste(outDir,SP,'_Timeseries_Site_',siteList[iSite],'.png',sep=''), width = 6, height = 2.5, units = "in", res=300)

      op <- par(mar=c(3, 6, 1, 1) + 0.1,mgp = c(2,1,0))
      xlabStr <- "Date"
    }else{
      png(paste(outDir,SP,'_Timeseries_Site_',siteList[iSite],'.png',sep=''), width = 6, height = 2,  units = "in", res = 300)

      op <- par(mar=c(1, 6, 1, 1) + 0.1,mgp = c(2,1,0))
      xlabStr <- ""
    }

    plot(AcOnlySegments$date[AcOnlySegments$siteNum==iSite],
         AcOnlySegments$Density[AcOnlySegments$siteNum==iSite],
         ylab = expression(atop(Estimated ~ Density,(animals/1000 ~ km^{2}))),
         xlab = xlabStr, pch = 1,col = alpha("black", 0.4),cex=.6,
         xlim = c(min(AcOnlySegments$date,na.rm = TRUE),max(AcOnlySegments$date,na.rm = TRUE)),
         ylim= c(0,max(AcOnlySegments$Density,na.rm=TRUE)))
    for (iRect in 1:length(deplGapsStart[[iSite]])){
        rect(deplGapsStart[[iSite]][iRect],0,deplGapsEnd[[iSite]][iRect],
             ceil(max(AcOnlySegments$Density, na.rm=TRUE)),
             border = "gray",col='gray')
    }
    points(AcOnlySegments$date[AcOnlySegments$siteNum==iSite],
         AcOnlySegments$Density[AcOnlySegments$siteNum==iSite],pch = 1,cex=.6,col = alpha("black", 0.4))

    text(x=max(AcOnlySegments$date,na.rm = TRUE),
         y = max(AcOnlySegments$Density*.95,na.rm = TRUE,
                 cex=1, adj = c(0,0)), labels = siteList[iSite])
    dev.off()
  }
}
