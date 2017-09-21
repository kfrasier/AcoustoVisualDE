library(raster)
library(rgdal)
library(lubridate)
library(geosphere)
library(sp)
# read in table of eddy centroid locations
eddyCentroidsTable <- read.csv('E:/NASData/Eddy/eddyCentroids.txt')

# Read in prediction template and change to lat/lon so that the distance 
# to centroid can be calculated for each cell location
predTemplate <- raster('E:/NASData/AcoustoVisualDE/Prediction_template/Predictiontemplate.tif')
crs(predTemplate)
r_pts <- rasterToPoints(predTemplate, spatial=TRUE) 
geo_projection <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
r_pts_tr <- spTransform(r_pts, CRS(geo_projection)) 
crs(r_pts_tr)

r_pts_tr@data <- data.frame(r_pts_tr@data, long = coordinates(r_pts_tr)[,1],
                         lat=coordinates(r_pts_tr)[,2])   
eddyCentroidsNegative <- eddyCentroidsTable[which(eddyCentroidsTable$cyclonic_type==-1),]
# Turn eddy dates into real dates,
eddy_dates <- as.POSIXct(strptime(eddyCentroidsNegative$obsdate,"%m/%d/%Y"),tz="GMT")

# for each month, find the mean min dist to centroid for each cell of prediction template
eddyLatLons <- eddyCentroidsNegative[c('longitude','latitude')]
eddyLatLons$longitude<- 180-eddyLatLons$longitude
yearList <- sort(unique(year(eddy_dates)))

minNegEddyDist <- matrix(NA, nrow(r_pts_tr),12)
for (iC in 1:nrow(r_pts_tr)){
  for (iM in 1:12){
    min_dist_set_year <- rep(NA, length(yearList))
    in_this_month <- which(month(eddy_dates) == iM)
    for (iY in 1:length(yearList)){
      this_year <- yearList[iY]
      in_this_year <- which(year(eddy_dates[in_this_month])==this_year)
      # for each day in month find nearest eddy to point, then average?
      day_set <- sort(unique(day(eddy_dates[in_this_month[in_this_year]])))
      
      # Initialize vector to store min eddy dist for each day
      min_dist_set_day <- rep(NA, length(day_set))
      for (iD in 1:length(day_set)){
        this_day <- day_set[iD]
        
        on_this_day <- which(day(eddy_dates[in_this_month[in_this_year]])== this_day)
        min_dist_set_day[iD]<- min(distm(r_pts_tr@data[iC,c(2,3)],eddyLatLons[in_this_month[in_this_year][on_this_day],]))
        
      }
      
      min_dist_set_year[iY]<- mean(min_dist_set_day) # mean across all days in month
      
    }
    minNegEddyDist[iC,iM] <- mean(min_dist_set_year) 
  }
  cat(paste0('Done with cell ', iC, ' of ', nrow(r_pts_tr),'\n'))
  
}

write.csv(minNegEddyDist,'E:/NASData/Eddy/minNegEddyDistTable.txt',na = "NA",row.names = FALSE,quote = FALSE)

janRaster <-predTemplate
janRaster[,]<-minNegEddyDist[,1]/1000
febRaster <-predTemplate
febRaster[,]<-minNegEddyDist[,2]/1000
marRaster <-predTemplate
marRaster[,]<-minNegEddyDist[,3]/1000
aprRaster <-predTemplate
aprRaster[,]<-minNegEddyDist[,4]/1000
mayRaster <-predTemplate
mayRaster[,]<-minNegEddyDist[,5]/1000
junRaster <-predTemplate
junRaster[,]<-minNegEddyDist[,6]/1000
julRaster <-predTemplate
julRaster[,]<-minNegEddyDist[,7]/1000
augRaster <-predTemplate
augRaster[,]<-minNegEddyDist[,8]/1000
sepRaster <-predTemplate
sepRaster[,]<-minNegEddyDist[,9]/1000
octRaster <-predTemplate
octRaster[,]<-minNegEddyDist[,10]/1000
novRaster <-predTemplate
novRaster[,]<-minNegEddyDist[,11]/1000
decRaster <-predTemplate
decRaster[,]<-minNegEddyDist[,12]/1000

writeRaster(janRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_janClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(febRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_febClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(marRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_marClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(aprRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_aprClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(mayRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_mayClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(junRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_junClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(julRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_julClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(augRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_augClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(sepRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_sepClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(octRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_octClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(novRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_novClimate.tif',
            format = 'GTiff',overwrite=TRUE)
writeRaster(decRaster,
            filename='E:/NASData/AcoustoVisualDE/EddyDist/NegEddyDist_decClimate.tif',
            format = 'GTiff',overwrite=TRUE)