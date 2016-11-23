# merge sighting sheets

visData1 <- read.csv("E:/NASData/Mammal_Sightings_9294.csv",na.strings=c(""," ","NA"))
visData2 <- read.csv("E:/NASData/Mammal_Sightings_9601.csv",na.strings=c(""," ","NA"))
visData3 <- read.csv("E:/NASData/GU0302_All_Sightings.csv",na.strings=c(""," ","NA"))
visData4 <- read.csv("E:/NASData/GU0402_All_Sightings.csv",na.strings=c(""," ","NA"))
visData5 <- read.csv("E:/NASData/GU0903_All_Sightings.csv",na.strings=c(""," ","NA"))
visData6 <- read.csv("E:/NASData/GU1202_Sightings.csv",na.strings=c(""," ","NA"))
visData7 <- read.csv("E:/NASData/GU1404_Sightings.csv",na.strings=c(""," ","NA"))


visData <-NULL
visData$date <-c(as.Date(visData1$date_,"%m/%d/%Y"), as.Date(visData2$date_,"%m/%d/%Y"),
                 as.Date(visData3$date_,"%m/%d/%Y"), as.Date(visData4$date_,"%m/%d/%Y"),
                 as.Date(visData5$date_,"%m/%d/%Y"), as.Date(visData6$date_,"%m/%d/%Y"),
                 as.Date(visData7$date_,"%m/%d/%Y"))

visData$sighting <-c(visData1$sighting,visData2$sighting,visData3$sighting,visData4$sighting,visData5$sighting,
                     visData6$sighting,visData7$sighting)

visData$commonname <- unlist(list(visData1$commonname,visData2$commonname,visData3$commonname,visData4$commonname,
                        visData5$commonname,visData6$commonname,visData7$commonname))

visData$spetype <- unlist(list(visData1$spetype,visData2$spetype,visData3$spetype,visData4$spetype,
                     visData5$spetype,visData6$spetype,visData7$spetype))

visData$size <- c(visData1$size,visData2$size,visData3$size,visData4$size,
                  visData5$size,visData6$size,visData7$size)

visData$effort <- c(visData1$effort,visData2$effort,visData3$effort,visData4$effort,
                   visData5$effort,visData6$effort,visData7$effort)

visData$calves <- c(visData1$calves,visData2$calves,visData3$calves,visData4$calves,
                     visData5$calves,visData6$calves,visData7$calves)

visData$dist <- c(visData1$dist,visData2$dist,visData3$dist,visData4$dist,
                  visData5$dist,visData6$dist,visData7$dist)

visData$distunit <- unlist(list(visData1$distunit,visData2$distunit,visData3$distunit,visData4$distunit,
                      visData5$distunit,visData6$distunit,visData7$distunit))

visData$relbear <- c(visData1$relbear,visData2$relbear,visData3$relbear,visData4$relbear,
                     visData5$relbear,visData6$relbear,visData7$relbear)

visData$beardir <- c(visData1$beardir,visData2$beardir,visData3$beardir,visData4$beardir,
                     visData5$beardir,visData6$beardir,visData7$beardir)

visData$dist_m <- c(visData1$dist_m,visData2$dist_m,visData3$dist_m,visData4$dist_m,
                    visData5$dist_m,visData6$dist_m,visData7$dist_m)

visData$boatlat <- c(visData1$boatlat,visData2$boatlat,visData3$boatlat,visData4$boatlat,
                     visData5$boatlat,visData6$boatlat,visData7$boatlat)

visData$boatlon <- c(visData1$boatlon,visData2$boatlon,visData3$boatlon,visData4$boatlon,
                     visData5$boatlon,visData6$boatlon,visData7$boatlon)

visData$transect_distm <- c(visData1$transect_distm,visData2$transect_distm,visData3$transect_distm,visData4$transect_distm,
                            visData5$transect_distm,visData6$transect_distm,visData7$transect_distm)

visData$transect <- c(visData1$transect,visData2$transect,visData3$transect,visData4$transect,
                      visData5$transect,visData6$transect,visData7$transect)

visData$seastate <- c(visData1$seastate,visData2$seastate,visData3$seastate,visData4$seastate,
                      visData5$seastate,visData6$seastate,visData7$seastate)

visData$vis <- c(visData1$vis,visData2$vis,visData3$vis,visData4$vis+1,
                 visData5$vis+1,visData6$vis+1,visData7$vis+1)

visData$weather <- c(visData1$weather,visData2$weather,visData3$weather,visData4$weather,
                     visData5$weather,visData6$weather,visData7$weather)

visData$glare <- c(visData1$glare,visData2$glare,visData3$glare,visData4$glare,
                   visData5$glare,visData6$glare,visData7$glare)


visData$cloud <- c(visData1$cloud,visData2$cloud,visData3$cloud,visData4$cloud,
                    visData5$cloud,visData6$cloud,visData7$cloud)

visData5$swell[which(visData5$swell<0)]<-NA
visData$swell <- c(visData1$swell,visData2$swell,visData3$swell,visData4$swell,
                   visData5$swell,visData6$swell,visData7$swell)

visData$windspeed <- c(visData1$windspeed,visData2$windspeed,visData3$windspeed,visData4$windspeed,
                       visData5$windspeed,visData6$windspeed,visData7$windspeed)

visData$flying_height <- c(visData1$flying_height,visData2$flying_height,visData3$flying_height,visData4$flying_height,
                           visData5$flying_height,visData6$flying_height,visData7$flying_height)

visData$ship <- unlist(list(visData1$ship,visData2$ship,visData3$ship,visData4$ship,
                  visData5$ship,visData6$ship,visData7$ship))
visData <- as.data.frame(visData)

cList <- colnames(visData)
png('Merged_Sighting_Verification.png', width = 800, height = 500)
par(mfrow = c(5, 5),mar = c(2,1,2,1))
for (itr in cList){
  myData <- visData[,itr]
  plot(myData[order(visData$date)],main = itr)
}
dev.off()

save(visData, file = "Sightings_merged.Rdata")