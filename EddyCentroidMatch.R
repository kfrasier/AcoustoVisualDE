library(geosphere)
load('E:/NASData/ModelData/Gg/setup_info_Gg.Rdata')
eddyCentroidsTable <- read.csv('E:/NASData/Eddy/eddyCentroids.txt')
eddyCentroidsTable$longitude<- 180-eddyCentroidsTable$longitude
eddyCentroidsTablePruned<-eddyCentroidsTable[which(eddyCentroidsTable$amplitude>=10),]
eddyDates <- as.POSIXct(strptime(eddyCentroidsTable$obsdate,"%m/%d/%Y"),tz="GMT")

acSegmentsAll <- read.csv(acousticSegFile, header = TRUE,na.strings=c(""," ","NA","-99999","-9999","NaN"))
acSegmentsAll$XLSDATE <- as.POSIXct(acSegmentsAll$XLSDATE,"%Y-%m-%d",tz = "GMT")
acSegmentsAll$EddyCenterDist <-NA
acSegmentsAll$PosEddyCenterDist <-NA
acSegmentsAll$NegEddyCenterDist <-NA
for (iR in 1:nrow(acSegmentsAll)){
  segDate <- acSegmentsAll$XLSDATE[iR]
  onThisDay <- which(segDate ==eddyDates)
  if (length(onThisDay)!=0){
    acSegmentsAll$EddyCenterDist[iR] <- min(distm(acSegmentsAll[iR,c('LONG','LAT')],
                                                  eddyCentroidsTable[onThisDay,c('longitude','latitude')]))
    posDir <- which(eddyCentroidsTable$cyclonic_type[onThisDay]==1)
    negDir <- which(eddyCentroidsTable$cyclonic_type[onThisDay]==-1)
    acSegmentsAll$PosEddyCenterDist[iR] <- min(distm(acSegmentsAll[iR,c('LONG','LAT')],
                                                     eddyCentroidsTable[onThisDay[posDir],c('longitude','latitude')]))
    acSegmentsAll$NegEddyCenterDist[iR] <- min(distm(acSegmentsAll[iR,c('LONG','LAT')],
                                                     eddyCentroidsTable[onThisDay[negDir],c('longitude','latitude')]))
    
  }
}

write.csv(acSegmentsAll,acousticSegFile,na = "NA",row.names = FALSE,quote = FALSE)


visSegments <- read.csv(visSegmentsFile, header = TRUE,sep = ",",na.strings=c(""," ","NA","-99999.0000","-99999"))
visSegments$date <- as.POSIXct(strptime(visSegments$date,"%Y-%m-%d"),tz="GMT")
visSegments$EddyCenterDist <- NA
visSegments$PosEddyCenterDist <- NA
visSegments$NegEddyCenterDist <- NA
for (iR in 1:nrow(visSegments)){
  visSegDate <- visSegments$date[iR]
  onThisDay <- which(visSegDate ==eddyDates)
  if (length(onThisDay)!=0){
    visSegments$EddyCenterDist[iR] <- min(distm(visSegments[iR,c('Long','Lat')],
                                                eddyCentroidsTable[onThisDay,c('longitude','latitude')]))
    posDir <- which(eddyCentroidsTable$cyclonic_type[onThisDay]==1)
    negDir <- which(eddyCentroidsTable$cyclonic_type[onThisDay]==-1)
    visSegments$PosEddyCenterDist[iR] <- min(distm(visSegments[iR,c('Long','Lat')],
                                                     eddyCentroidsTable[onThisDay[posDir],c('longitude','latitude')]))
    visSegments$NegEddyCenterDist[iR] <- min(distm(visSegments[iR,c('Long','Lat')],
                                                     eddyCentroidsTable[onThisDay[negDir],c('longitude','latitude')]))
  }
}
write.csv(visSegments,visSegmentsFile,na = "NA",row.names = FALSE,quote = FALSE)


