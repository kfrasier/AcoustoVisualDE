#' Function to plot acoustic timeseries data and save to image. 
#' @export 
#' @examples 
#' plot.timeseries() 

plot.timeseries <-function(siteList,outDir,AcOnlySegments){
  
  for (iSite in 1:length(siteList)){
    
    if (iSite == length(siteList)){
      png(paste(outDir,SP,'_Timeseries_Site_',siteList[iSite],'.png',sep=''), width = 480, height = 130+30)
      
      op <- par(mar=c(3, 6, 0, 1) + 0.1,mgp = c(2,1,0))
      xlabStr <- "Date"
    }else{
      png(paste(outDir,SP,'_Timeseries_Site_',siteList[iSite],'.png',sep=''), width = 480, height = 130)
      
      op <- par(mar=c(1, 6, 0, 1) + 0.1,mgp = c(2,1,0))
      xlabStr <- ""
    }
    
    plot(AcOnlySegments$date[AcOnlySegments$siteNum==iSite],
         AcOnlySegments$Density[AcOnlySegments$siteNum==iSite],
         ylab = expression(atop(Estimated ~ Density,(animals/1000 ~ km^{2}))),
         xlab = xlabStr, 
         xlim = c(min(AcOnlySegments$date,na.rm = TRUE),max(AcOnlySegments$date,na.rm = TRUE)))
    # rect()
    text(x=max(AcOnlySegments$date,na.rm = TRUE),
         y = max(AcOnlySegments$Density[AcOnlySegments$siteNum==iSite]*.95,na.rm = TRUE,
                 cex=1.2, adj = c(0,0)), labels = siteList[iSite])
    dev.off()
  }
}