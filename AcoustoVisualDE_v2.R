# Acoustovisual density estimation
library(mrds)
library(lubridate)
library(magic)
library(mgcv)
library(openair)
library(prodlim)
library(psych)
library(pracma)
library(plotrix)
library(HabitatProject)
library(rgdal)
library(raster)

## Read set up file for species of choice.
# NOTE: if you have changed the setup info, re-run "setup_info_[your species here].R" before running this
load('E:/NASData/ModelData/Zc/setup_info_Zc.Rdata')

# Set up directories
outDir <- file.path("E:/NASData/ModelData",SP,"/")
setwd(outDir)


############### Begin Calculations ################
graphics.off()
closeAllConnections()

############### Load & Prune Acoustic Data ################
cat("Loading acoustic data \n")
if (matchACSegs){
  # Load acoustic segments and densities
  acSegmentsAll <- read.csv(acousticSegFile, header = TRUE,na.strings=c(""," ","NA","-99999","-9999","NaN"))
  acDensityAll <- read.csv(acousticDensityFile, header = TRUE,na.strings=c(""," ","NA","-99999","-9999","NaN"))
  acSegmentsAll$XLSDATE <- as.POSIXct(acSegmentsAll$XLSDATE,"%Y-%m-%d",tz = "GMT")
  acDensityAll$xlsDate <- as.POSIXct(strptime(acDensityAll$xlsDate,"%m/%d/%Y"),tz = "GMT")#"%m/%d/%Y %H:%M"
  nCol <- length(colnames(acDensityAll))
  
  keepPoints <- which(acDensityAll$xlsDate >= "2011-01-01" & acDensityAll$xlsDate < "2014-01-01")

  acDensityAll <- acDensityAll[keepPoints,]
  
  acSegLatFloor <- floor(acSegmentsAll$LAT *1000)/1000
  acDensLatFloor <- floor(acDensityAll$Lat *1000)/1000
  
  # Match segments to density datapoints
  covarNames <- names(acSegmentsAll[5:length(names(acSegmentsAll))])
  acDensityAll[,covarNames] <- NA
  for (iR in 1:nrow(acDensityAll)){
    # find all the segments with matching latitudes
    
    goodLat <- which(acSegLatFloor %in% acDensLatFloor[iR])
    densDate <- acDensityAll$xlsDate[iR]
    bestMatch <- goodLat[which.min(abs(densDate-acSegmentsAll$XLSDATE[goodLat]))]
    
    if (length(bestMatch)!=0){
      rowMatch <- acSegmentsAll[bestMatch,5:(length(names(acSegmentsAll)))]
      acDensityAll[iR,((nCol+1):(length(rowMatch)+nCol))] <- rowMatch
    }else{
      acDensityAll[iR,((nCol+1):(length(covarNames)+nCol))] <- NA
    }
    
    if (iR %% 1000 == 0){
      cat(paste0("done with entry ", iR , " of ", nrow(acDensityAll), "\n", collapse = ""))
    }
  }
  save(acDensityAll, file = acDensityFile)
  rm(acSegmentsAll)

}else {
  load(acDensityFile)
}

# Exclude partial weeks, and extract the right density estimate type (cue or group)
#fullWeeks <- which(acSegmentsAll$PartialWeek == 0)
#acSegmentsFull <- acSegmentsAll#[fullWeeks,]

####################### Load & Prune Visual Data #########################
cat("Loading visual data\n")

load(visDataFile)
# Truncate visual data to exclude sightings outside the 200m contour polygon
visDataDF <- data.frame(visData)
pred_polygon <- readOGR('E:/NASData/AcoustoVisualDE/Prediction_template/prediction_polygon.shp')
coordinates(visDataDF) <- ~boatlon + boatlat
crs(visDataDF) <- crs("+proj=longlat +datum=WGS84")
pred_polygon_proj <- spTransform(pred_polygon,crs(visDataDF))
visDataDF_proj_crop <- visDataDF[pred_polygon_proj,]
visDataDF_crop <- as.data.frame(visDataDF_proj_crop)
# visDataDF_crop<-visDataDF

visSegments <- read.csv(visSegmentsFile, header = TRUE,na.strings=c(""," ","NA","-99999.0000","-99999"))
visSegments$date <- as.POSIXct(strptime(visSegments$date,"%Y-%m-%d"),tz="GMT")
## Process Visual data to determine detection probabilities and strip widths
visSpIdx <- NULL
for (spname in SPC_vis){
  visSpIdx <- c(visSpIdx,which(grepl(spname, visDataDF_crop$commonname)))
}

# prune out off-effort sightings
spIdxON <- visSpIdx[which(as.logical(visDataDF_crop$effort[visSpIdx]))]

# prune out sightings with an angle over 90 deg
spIdxON2 <- spIdxON[which(visDataDF_crop$relbear[spIdxON]<90)]

# Populate truncated column of on effort sightings of species of interest with zeros
visDataDF_crop$Truncated <- 1
visDataDF_crop$Truncated[spIdxON2] <- 0

####################### Fit Visual Survey Detection Functions #########################

if (runDetFuns){
  # for each visual platform
  nPlatform <- 1  # Initialize to start with first platform
  
  tDist <-NULL # store truncation distances
  bestModel <-NULL # store best model covariates
  bestKey <-NULL # store best model key
  detFunByPlatform <- NULL
  
  cat("Begin model fitting for visual data \n")
  for (i in PLC){
    ddfData <- NULL
    cat(paste("Fitting platform ", i,"\n"))
    
    # identify sightings assoicated with the a certain platform
    PLCspIdxOn <- spIdxON2[which(grepl(i,visDataDF_crop$ship[spIdxON]))] 
    
    # Make dataframe with the inputs that the ddf distance function wants
    ddfData$observer <- rep(1,length(PLCspIdxOn))
    ddfData$detected <- rep(1,length(PLCspIdxOn))
    ddfData$object <- (1:length(PLCspIdxOn))
    ddfData$distance <- visDataDF_crop$transect_distm[PLCspIdxOn]
    ddfData$size <- visDataDF_crop$size[PLCspIdxOn]
    ddfData$seastate <- visDataDF_crop$seastate[PLCspIdxOn]
    ddfData$swell <- visDataDF_crop$swell[PLCspIdxOn]
    ddfData$vis <- visDataDF_crop$vis[PLCspIdxOn]
    ddfData <- as.data.frame(ddfData)
    
    # Compute untruncated detection function
    cat("Calculating basic fit with non-truncated data and half-normal key, no covariates.\n")
    
    detFun_noTrunc <-ddf(method = 'ds', dsmodel =~ mcds(key = 'hn', formula = ~ 1),
                         data = ddfData, meta.data = list(binned=F,left=0),
                         control = list(optimx.maxit = 20))
    # detFun_noTrunc2 <-ds(as.data.frame(ddfData),key='hn')
    
    
    # Make output plots
    cat("Saving plots\n")
    png(paste(outDir, SP,'sightnoTrunc_',i,'.png',sep=''), width = 800, height = 500)
    par(mfrow=c(1,2))
    plot(detFun_noTrunc)
    qqplot.ddf(detFun_noTrunc,plot=TRUE)
    dev.off()
    
    
    # Compute truncation distance by removing highest 5% of distances
    tDist[nPlatform] <- quantile(ddfData$distance,.95,na.rm = TRUE)
    cat(paste("Truncation distance for platform ", i, "=",  round(tDist[nPlatform],2), "m \n"))
    
    
    # Iterate over detection functions with various adjustments and orders, and identify AIC for each
    # list of key funs to try: 
    keyListInit = c(  'hn',   'hr',  'hr','unif')
    adjInit      = c('none','none','poly','none')
    adjOrderInit = c(     0,     0,     2,     0)
    detFun1 <- NULL
    aicList1 <- NULL
    keyList1 <- NULL
    adjStr <- NULL
    cat("Fitting detection functions with adjustments \n")
    
    dI <- 1
    for (i1 in 1:length(keyListInit)){
      df <-NULL
      if (grepl('none',adjInit[[i1]])){
        df <- ddf(method ='ds', dsmodel =~ mcds(key = keyListInit[[i1]], formula = ~ 1),
                  data = ddfData, meta.data = list(binned=F, width=tDist[nPlatform], left=0))
        #         df <- ds(ddfData, truncation = tDist[nPlatform], order = NULL, transect = "line", key = keyListInit[[i1]],
        #                   monotonicity = "weak")
        
      } else {
        df <-ddf(method ='ds', dsmodel =~ mcds(key = keyListInit[[i1]], formula = ~ 1,
                                               adj.series = adjInit[[i1]], adj.order = adjOrderInit[[i1]]), data = as.data.frame(ddfData),
                 meta.data = list(binned=F, width=tDist[nPlatform],left=0))
        #         df <- ds(ddfData, truncation = tDist[nPlatform], transect = "line", key = keyListInit[[i1]],
        #                  adjustment = adjInit[i1], order = adjOrderInit[i1],
        #                  monotonicity = "weak")
      }
      
      
      if (is.null(df)){
        cat(paste0("Model did not converge: Key = ",iKey, "; key = ",keyList1[iI],"\n", collapse = ""))
        
      }else {
        detFun1[[dI]] <-df
        
        aicList1[dI] <- df$criterion
        keyList1[dI] <- keyListInit[i1]
        adjStr[dI] <- paste(adjInit[i1],adjOrderInit[i1])
        
        cat(paste0("Model result ", i1,"; Key = ",keyList1[i1],
                   "; Adjustment = ",adjInit[i1], "; Order = ",adjOrderInit[i1],"\n", collapse = ""))
        cat(paste0("AIC = ",  round(aicList1[i1], digits=2),"\n", collapse = ""))
        dI <- dI+1
      }
    }
    
    
    # #### Fitting detection functions with covariates
    # # Iterate over detection functions with covariates and identify AIC for each
    # covarList <- c("size", "seastate",'vis', 'swell')
    # # list of keys: hn, hr
    # keyCovar= c('hn', 'hr')
    # CI <- 1
    # aicList2 <- NULL # store AIC scores
    # keyList2 <- NULL # store keys scores
    # cSetStr <- NULL # store the covariate formulas
    # detFun2<- NULL
    # adjList2 <-NULL
    # 
    # # iterate over the key options
    # for (iKey in keyCovar){
    #   # iterate over covariate combinations - using "combn" to come up with different the combinations of covariates
    #   for (i2 in 1 : length(covarList)){
    #     covarSet <- combn(covarList,i2)
    #     nSets <- dim(covarSet)
    #     # for each set of covariates, fit a model
    #     for (i3 in 1:nSets[2]){
    #       cSet <- covarSet[,i3]
    #       # nPlus <- nSets[1]-1
    #       cSetStr[CI] <- paste('~', paste0(cSet, collapse = " + "))
    #       # sometimes models do not converge, use try() to avoid crash if a model fails
    #       dF <- NULL
    #       try(dF <- ddf(method='ds',dsmodel=~mcds(key=iKey,  formula = cSetStr[CI]),
    #                     data = ddfData,
    #                     meta.data = list(binned=F, width=tDist[nPlatform],left=0)))
    #       
    #       if (is.null(dF)){
    #         cat(paste0("Model did not converge: Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
    #         
    #       }else {
    #         detFun2[[CI]] <-dF
    #         adjList2[CI]<- cSetStr[CI]
    #         aicList2[CI] <- detFun2[[CI]]$criterion
    #         keyList2[CI] <- iKey
    #         cat(paste0("Model result ", CI+i1, ": Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
    #         cat(paste("AIC =",  round(aicList2[CI], digits=2),"\n", collapse = ""))
    #         CI <- CI+1
    #       }
    #     }
    #   }
    # }
    # 
    # cat("Done fitting models")
    
    # Put all combinations together, and see which one has the lowest AIC
    aicList<-c(aicList1)#,aicList2) 
    keyList <- c(keyList1)#,keyList2)
    adjList <- c(adjStr)#,adjList2)
    detFun <- c(detFun1)#,detFun2)
    
    ddfOut <- data.frame(model = adjList, key = keyList, aic = aicList)
    
    # Best model is...
    bestModelIdx <- which(ddfOut$aic == min(ddfOut$aic, na.rm = TRUE), arr.ind = TRUE)
    bestModel[[nPlatform]] <- adjList[bestModelIdx]
    bestKey[[nPlatform]] <- keyList[bestModelIdx]
    cat(paste("Best model for Platform ", i,":\n"))
    cat(paste("Key = ", bestKey[[nPlatform]], "; Adjustment =",  bestModel[[nPlatform]],"\n"))
    
    
    # Make output plot of best model
    cat("Saving plots and summaries \n")
    png(paste(outDir, SP,'sightwTrunc_',i,'.png',sep=''), width = 800, height = 500,pointsize = 16)
    par(mfrow=c(1,2))
    plot(detFun[[bestModelIdx]], main = paste('model = ',bestModel[[nPlatform]], '; key = ', bestKey[[nPlatform]]))
    qqplot.ddf(detFun[[bestModelIdx]],plot=TRUE)
    dev.off()
    
    # Output summary text to txt file
    sink(paste(outDir, SP,'sightwTrunc_',i,'.txt',sep=''))
    print(summary(detFun[[bestModelIdx]]))
    print(ddf.gof(detFun[[bestModelIdx]]))
    sink()
    
    ## save best model
    save(detFun, bestModelIdx, file = paste(SP,'sightwTrunc_',i,'.Rdata',sep=''))
    detFunByPlatform[[nPlatform]] <- detFun[[bestModelIdx]]
    
    
    nPlatform <- nPlatform + 1
    cat(paste("Done fitting models for platform ", i,"\n"))
  }
  
  cat("Done fitting models \n")
  save(tDist, detFunByPlatform, file = detFunFile)
  
} else{
  load(detFunFile)
}

####################### Match Sightings with Transects #########################

# Put effective strip widths calculated above into Segments table
visSegments["ESW"] <- 0
for (iP in 1:length(PLC)){
  thisSegList <-which(visSegments$ship==PLC[iP])
  visSegments$ESW[thisSegList] <- predict(detFunByPlatform[[iP]],
                                          esw=TRUE)$fitted[1]
}

# Get rid of off effort segments
visSeg_OnEffort <- visSegments[which(visSegments$effort==1),]

# Tally encounters by segment
prunedSightings <- visDataDF_crop[which(visDataDF_crop$Truncated==0),] # get all of the non-truncated sightings

# Assign segment to each sighting
for (iSight in 1:length(prunedSightings$date)){
  sightDate <- as.POSIXct(strptime(prunedSightings$date[iSight],"%Y-%m-%d"),tz="GMT")
  # cat(paste0('Original date = ', prunedSightings$date[iSight],'\n'))
  # cat(paste0('POSIXct date = ', sightDate,'\n'))
  onThisDay <- which(visSeg_OnEffort$date == sightDate)
  # cat(paste0('matching day = ', onThisDay,'\n'))
  if (length(onThisDay)>0) {
    minIdx <- which.min(rowSums((visSeg_OnEffort[onThisDay,c('Lat','Long')]-
              matrix(as.numeric(rep(prunedSightings[iSight,c('boatlat','boatlon')],each=length(onThisDay)),ncol=2)))^2))
    prunedSightings$Segment[iSight] <- onThisDay[minIdx]
    
  }else { # handle case where there is no match (why would this happen?)
    cat(paste("Warning: Missing effort segment for sighting on", sightDate,"\n"))
    prunedSightings$Segment[iSight] <- NaN
  }
}

#adjust encounters for G0
# for (i in PLC){
#   thisSet <- which(visSeg_OnEffort$ship = i)
#   visSeg_OnEffort$SpEncounter_G0adj[thisSet] <- visSeg_OnEffort$sp_count[thisSet]
# }

cat("Associating sightings with transect segments\n")

# put that info into the segments table
visSeg_OnEffort$sp_count <- 0  # will hold number of animals
visSeg_OnEffort$sp_present <- 0 # will hold 1/0 for presence absence

# for each segment that had a sighting
for (iSeg in 1:length(prunedSightings$Segment)){
  
  # Add the number of animals observed, 
  # this is cumulative in case multiple sightings occured on one transect
  visSeg_OnEffort$sp_count[prunedSightings$Segment[iSeg]] <- 
    visSeg_OnEffort$sp_count[prunedSightings$Segment[iSeg]] + 
    prunedSightings$size[iSeg]

  # also set presence equal to 1
  visSeg_OnEffort$sp_present[prunedSightings$Segment[iSeg]] <- 1
}

# account for G0 in encounters
visSeg_OnEffort$sp_count_g0adj <- visSeg_OnEffort$sp_count/(visG0* detFun[[bestModelIdx]]$fitted[1])

# Estimate surveyed area
visSeg_OnEffort$EffectiveArea <- (2*tDist/1000)*(visSeg_OnEffort$SegmentLength/1000)
visSeg_OnEffort$Density <- (visSeg_OnEffort$sp_count_g0adj)/visSeg_OnEffort$EffectiveArea


############## Form Acoustic and Visual Covariate Dataframes #####################

AcOnlySegments <- NULL
yearListTemp <- as.numeric(format(acDensityAll$xlsDate,"%Y"))
siteYear <- NULL
nAc <- length(acDensityAll$HYCOM_MAG_100)
siteYear$Year <-yearListTemp
siteYear$Site <-acDensityAll$Site
siteYear <- as.data.frame(siteYear)
uSiteYear <- unique((siteYear))
# AcOnlySegments$Density <- rep(NA,times = nAc)
# for (uR in 1:nrow(uSiteYear)){
#   # Normalize acoustic density estimated by deployment
#   thisSet <- which(as.logical(row.match(siteYear,uSiteYear[uR,])))
#   thisSet_gt0 <- which(acDensityAll$meanDensity[thisSet]>0)
#   quant95 <-quantile(acDensityAll$meanDensity[thisSet[thisSet_gt0]],probs = .95,na.rm = TRUE)
AcOnlySegments$Density <- acDensityAll$meanDensity/(pi*r_sp^2)*1000#[thisSet]/quant95

AcOnlySegments$date <- acDensityAll$xlsDate #date
AcOnlySegments$Numeric_date <- (as.numeric(acDensityAll$xlsDate)-min(as.numeric(acDensityAll$xlsDate)))/100
AcOnlySegments$lat <- acDensityAll$Lat
AcOnlySegments$long <- acDensityAll$Long

AcOnlySegments$Category<- rep(2,length(acDensityAll$Long))
AcOnlySegments$SST <- acDensityAll$SST_DAILY_CMC
AcOnlySegments$SSH <- acDensityAll$SSH_DAILY_AVISO
AcOnlySegments$CHL <- acDensityAll$CHL_8DAY_NASA
#AcOnlySegments$HYCOM_QTOT <- acDensityAll$HYCOM_QTOT
AcOnlySegments$HYCOM_MLD <- acDensityAll$HYCOM_MLD
AcOnlySegments$HYCOM_EMP <- acDensityAll$HYCOM_EMP
AcOnlySegments$HYCOM_DIR_0 <- acDensityAll$HYCOM_DIR_0
AcOnlySegments$HYCOM_DIR_100 <- acDensityAll$HYCOM_DIR_100
AcOnlySegments$HYCOM_SALIN_0 <- acDensityAll$HYCOM_SALINITY_0
AcOnlySegments$HYCOM_SALIN_100 <- acDensityAll$HYCOM_SALIN_100
AcOnlySegments$HYCOM_MAG_0 <- acDensityAll$HYCOM_MAG_0
AcOnlySegments$HYCOM_MAG_100 <- acDensityAll$HYCOM_MAG_100
AcOnlySegments$HYCOM_UPVEL_100 <- acDensityAll$HYCOM_UPVEL_100
AcOnlySegments$HYCOM_UPVEL_50 <- acDensityAll$HYCOM_UPVEL_50
AcOnlySegments$FrontDist_Cayula <- acDensityAll$FRONTDIST_CAYULA
AcOnlySegments$EddyDist <- acDensityAll$EddyDist
AcOnlySegments$Neg_EddyDist <- acDensityAll$NegEddyDist
AcOnlySegments$Pos_EddyDist <- acDensityAll$PosEddyDist

AcOnlySegments$Type <- rep(2,times = nAc)
AcOnlySegments$DayOfYear <- as.numeric(strftime(acDensityAll$xlsDate,"%j")) # day of year
AcOnlySegments$Weight <- 24*(AcTruncDist^2)*pi # 24 hours, 5km radius
AcOnlySegments$Site<- acDensityAll$Site
AcOnlySegments <- as.data.frame(AcOnlySegments)
# Make vector indicating deployment categories based on lat/long
myLatLon = acDensityAll[c(4,5)]
# # uLatLon <- unique(myLatLon)
uLatLon.site <- unique(floor(myLatLon)) # indentify distinct sites based on lat/long
# uLatLon.deployment <- unique(floor(myLatLon*10000))/10000 # indentify distinct deployments based on lat/long
nRows <- length(acDensityAll[,4])
AcOnlySegments$siteNum <- rep(NA,times = nRows) # will hold site label
# AcOnlySegments$fac2 <- rep(NA,times = nRows) # will hold deployment label
uLatLon.site <- uLatLon.site [!is.na(uLatLon.site[,1]),]
for (uR in 1:nrow(uLatLon.site)){
  thisSet1 <- which(as.logical(row.match(floor(myLatLon),uLatLon.site[uR,])))
  AcOnlySegments$siteNum[thisSet1] <- uR
## indentify distinct deployments at site based on lat/long
#   uLatLon.deployment <- unique(floor(myLatLon[thisSet1,]*10000))/10000 
#   
#   for (uV in 1:nrow(uLatLon.deployment)){
#     thisSet2 <- which(as.logical(row.match(floor(myLatLon[thisSet1,]*10000)/10000,uLatLon.deployment[uV,])))
#     AcOnlySegments$fac2[thisSet1[thisSet2]] <-uV
#   }
}

VisOnlySegments <- NULL
VisOnlySegments$date <- as.POSIXct(strptime(visSeg_OnEffort$date_Converted,"%m/%d/%Y"),tz="GMT")#date
VisOnlySegments$Numeric_date <- (as.numeric(VisOnlySegments$date)-min(as.numeric(VisOnlySegments$date)))/100
VisOnlySegments$lat <- visSeg_OnEffort$Lat
VisOnlySegments$long <- visSeg_OnEffort$Long
VisOnlySegments$Category<- rep(1,length(visSeg_OnEffort$Long))
# mergedSegments$ESW <- c(acDensityAll$BW_ESW)
VisOnlySegments$Density <- visSeg_OnEffort$Density*1000
  #quantile(visSeg_OnEffort$Density[which(visSeg_OnEffort$Density>0)],probs = .95,na.rm = TRUE)
VisOnlySegments$SST <- visSeg_OnEffort$SST_daily_CMC_L4_GLOB
VisOnlySegments$SSH <- visSeg_OnEffort$SSH_daily_aviso_double
VisOnlySegments$CHL <- visSeg_OnEffort$CHl_8Day_NASA
# VisOnlySegments$HYCOM_QTOT <- visSeg_OnEffort$HYCOM_qTot
VisOnlySegments$HYCOM_MLD <- visSeg_OnEffort$HYCOM_mld
VisOnlySegments$HYCOM_EMP <- visSeg_OnEffort$HYCOM_emp
VisOnlySegments$HYCOM_DIR_0 <- visSeg_OnEffort$HYCOM_dir_0
VisOnlySegments$HYCOM_DIR_100 <- visSeg_OnEffort$HYCOM_dir_100
VisOnlySegments$HYCOM_SALIN_0 <- visSeg_OnEffort$HYCOM_salinity_0
VisOnlySegments$HYCOM_SALIN_100 <- visSeg_OnEffort$HYCOM_salin_100
VisOnlySegments$HYCOM_MAG_0 <- visSeg_OnEffort$HYCOM_mag_0
VisOnlySegments$HYCOM_MAG_100 <- visSeg_OnEffort$HYCOM_mag_100
VisOnlySegments$HYCOM_UPVEL_100 <- visSeg_OnEffort$HYCOM_upVel_100
VisOnlySegments$HYCOM_UPVEL_50 <- visSeg_OnEffort$HYCOM_UPVEL_50
VisOnlySegments$FrontDist_Cayula <- visSeg_OnEffort$FrontDist_Cayula
VisOnlySegments$EddyDist <- visSeg_OnEffort$EddyDist
VisOnlySegments$Neg_EddyDist <- visSeg_OnEffort$NegEddyDist
VisOnlySegments$Pos_EddyDist <- visSeg_OnEffort$PosEddyDist

nVis <- length(visSeg_OnEffort$HYCOM_upVel_100)
VisOnlySegments$Type <- rep(1,times = nVis)
VisOnlySegments$DayOfYear <- as.numeric(strftime(VisOnlySegments$date,format="%j")) # day of year
VisOnlySegments <- as.data.frame(VisOnlySegments)
visOrder <- order(VisOnlySegments$date)
VisOnlySegments <-VisOnlySegments[visOrder,] # make sure the segments are sequential 
# in case it matters for correlation structure.
VisOnlySegments$Weight <- (((visSeg_OnEffort$SegmentLength/1000)/18.52)*
                             visSeg_OnEffort$EffectiveArea)
      # transect length, divided by speed (10knots/hr = 18.52 km/hr), times estimated strip width, times length, *2

##################### Merge Visual and Acoustic Segments #####################
cat("Merging Visual and Acoustic Segments\n")

# Merge visual and acoustic segments into one big dataframe
mergedSegments <- NULL
mergedSegments$date <- c(VisOnlySegments$date,AcOnlySegments$date) #date
mergedSegments$Numeric_date <- (as.numeric(mergedSegments$date)-
                                  min(as.numeric(mergedSegments$date)))/100 #date
mergedSegments$lat <- c(visSeg_OnEffort$Lat,acDensityAll$Lat)
mergedSegments$long <- c(visSeg_OnEffort$Long,acDensityAll$Long)
mergedSegments$Category<- c(rep(1,length(visSeg_OnEffort$Long)),rep(2,length(acDensityAll$Long)))
# mergedSegments$ESW <- c(acDensityAll$BW_ESW)
mergedSegments$Density <- c(VisOnlySegments$Density,
                            AcOnlySegments$Density)
mergedSegments$SST <- c(visSeg_OnEffort$SST_daily_CMC_L4_GLOB,acDensityAll$SST_DAILY_CMC)
mergedSegments$SSH <- c(visSeg_OnEffort$SSH_daily_aviso_double,acDensityAll$SSH_DAILY_AVISO)
mergedSegments$CHL <- c(visSeg_OnEffort$CHl_8Day_NASA, acDensityAll$CHL_8DAY_NASA)
# mergedSegments$HYCOM_QTOT <- c(visSeg_OnEffort$HYCOM_qTot, acDensityAll$HYCOM_QTOT)
mergedSegments$HYCOM_MLD <- c(visSeg_OnEffort$HYCOM_mld, acDensityAll$HYCOM_MLD)
mergedSegments$HYCOM_EMP <- c(visSeg_OnEffort$HYCOM_emp, acDensityAll$HYCOM_EMP)
mergedSegments$HYCOM_DIR_0 <- c(visSeg_OnEffort$HYCOM_dir_0,acDensityAll$HYCOM_DIR_0)
mergedSegments$HYCOM_DIR_100 <- c(visSeg_OnEffort$HYCOM_dir_100,acDensityAll$HYCOM_DIR_100)
mergedSegments$HYCOM_SALIN_0 <- c(visSeg_OnEffort$HYCOM_salinity_0,acDensityAll$HYCOM_SALINITY_0)
mergedSegments$HYCOM_SALIN_100 <- c(visSeg_OnEffort$HYCOM_salin_100,acDensityAll$HYCOM_SALIN_100)
mergedSegments$HYCOM_MAG_0 <- c(visSeg_OnEffort$HYCOM_mag_0,acDensityAll$HYCOM_MAG_0)
mergedSegments$HYCOM_MAG_100 <- c(visSeg_OnEffort$HYCOM_mag_100,acDensityAll$HYCOM_MAG_100)
mergedSegments$HYCOM_UPVEL_100 <- c(visSeg_OnEffort$HYCOM_upVel_100,acDensityAll$HYCOM_UPVEL_100)
mergedSegments$HYCOM_UPVEL_50 <- c(visSeg_OnEffort$HYCOM_UPVEL_50,acDensityAll$HYCOM_UPVEL_50)
mergedSegments$FrontDist_Cayula <- c(visSeg_OnEffort$FrontDist_Cayula,acDensityAll$FRONTDIST_CAYULA)
mergedSegments$EddyDist <- c(visSeg_OnEffort$EddyDist,acDensityAll$EddyDist)
mergedSegments$Neg_EddyDist <- c(visSeg_OnEffort$NegEddyDist,acDensityAll$NegEddyDist)
mergedSegments$Pos_EddyDist <- c(visSeg_OnEffort$PosEddyDist,acDensityAll$PosEddyDist)

mergedSegments$Type <- c(rep(1,times = nVis),rep(2,times = nAc))
mergedSegments$DayOfYear <- as.numeric(strftime(mergedSegments$date,"%j")) # day of year
mergedSegments$Weight <- c(VisOnlySegments$Weight,AcOnlySegments$Weight)
mergedSegments <- as.data.frame(mergedSegments)


############ Calculate a covariance factor based on deployment (Acoustic) or Year (Visual) #############
# Trying two grouping options: fac1 just groups acoustic monitoring by site, this seems ideal but results in small sample size
# So fac2 groups by deployment, which seems less defensible, but might help models converge.
# Visuals are grouped by cruise in both cases.
myLatLon <- data.frame(AcOnlySegments$lat,AcOnlySegments$long)
uLatLon <- unique(myLatLon)
notNA <- which(!is.na(uLatLon[,1]))
uLatLon <- uLatLon[notNA,]
nRows <- length(AcOnlySegments[,1])
fac2 <- rep(NA,times = nRows)
for (uR in 1:nrow(uLatLon)){
  thisSet <- which(as.logical(row.match(myLatLon,(uLatLon[uR,]))))
  fac2[thisSet] <-uR
}
AcOnlySegments$fac2 <- fac2
AcOnlySegments$fac1 <- AcOnlySegments$siteNum


myYear <- as.numeric(strftime(VisOnlySegments$date ,"%Y"))
uYear <- unique(myYear)
nRows <- length(myYear)
VisOnlySegments$fac1 <- rep(NA,times = nRows)
for (uR in 1:length(uYear)){
  thisSet <- which(as.logical(match(myYear,uYear[uR])))
  VisOnlySegments$fac1[thisSet] <- uR + max(AcOnlySegments$fac1,na.rm = TRUE)
}
VisOnlySegments$fac2 <- VisOnlySegments$fac1

mergedSegments$fac1 <- as.numeric(c(VisOnlySegments$fac1,AcOnlySegments$fac1))
mergedSegments$fac2 <- as.numeric(c(VisOnlySegments$fac2,AcOnlySegments$fac2))

save(mergedSegments,VisOnlySegments,AcOnlySegments,tDist,
     file = paste(outDir,SP,'MergedData.Rdata',sep=''))

##################### Data Exploration, Transformation and Plotting #####################

# Explore the data, graphical output
# cat("Begin exploratory plot generation\n")  
# plot.timeseries(siteList,outDir,AcOnlySegments)
# 
# # histograms of missing data
# percFilled <- plot.missingdata(mergedSegments,covarList,paste0('AcousticAndVisual_',SP)) 
# percFilled <- plot.missingdata(AcOnlySegments,covarList,paste0('AcousticOnly_',SP)) 
# percFilled <- plot.missingdata(VisOnlySegments,covarList,paste0('VisualOnly_',SP)) 
# 
# Identify and clear problematic outliers
# outlierList <-which(mergedSegments$CHL< -10)
# mergedSegments$CHL[outlierList] <- NaN 
# outlierList <-which(mergedSegments$FrontDist_Cayula>800000)
# mergedSegments$FrontDist_Cayula[outlierList] <- NaN 
# outlierList <-which(mergedSegments$Density>10000)
# mergedSegments$Density[outlierList] <- NaN 
# 
# outlierList <-which(AcOnlySegments$CHL< -10)
# AcOnlySegments$CHL[outlierList] <- NaN 
# outlierList <-which(AcOnlySegments$FrontDist_Cayula>800000)
# AcOnlySegments$FrontDist_Cayula[outlierList] <- NaN 
# outlierList <-which(AcOnlySegments$Density>10000)
# AcOnlySegments$Density[outlierList] <- NaN 
# 
# outlierList <-which(VisOnlySegments$CHL<  -10)
# VisOnlySegments$CHL[outlierList] <- NaN 
# outlierList <-which(VisOnlySegments$FrontDist_Cayula>800000)
# VisOnlySegments$FrontDist_Cayula[outlierList] <- NaN 
# outlierList <-which(VisOnlySegments$Density>10000)
# VisOnlySegments$Density[outlierList] <- NaN 

# If you decide from the missing data plots that you want to restrict years going forward:
# yearListIdx = as.numeric(format(mergedSegments$date,"%Y"))
# yearListIdx_AcOnly = as.numeric(format(AcOnlySegments$date,"%Y"))
# yearListIdx_VisOnly = as.numeric(format(VisOnlySegments$date,"%Y"))
# 
# keepDates.train <- which(yearListIdx != 2009 & yearListIdx >= 2003 & yearListIdx <= 2012)
# keepDates.test <- which(yearListIdx == 2009 | yearListIdx == 2013)
# keepDates_AcOnly.train <- which(yearListIdx_AcOnly != 2009 & yearListIdx_AcOnly >= 2003 & yearListIdx_AcOnly <= 2012)
# keepDates_AcOnly.test <- which(yearListIdx_AcOnly == 2009 | yearListIdx_AcOnly == 2013)
# keepDatesVisOnly.train <- which(yearListIdx_VisOnly != 2009 & yearListIdx_VisOnly >= 2003 & yearListIdx_VisOnly <= 2012)
# keepDatesVisOnly.test <- which(yearListIdx_VisOnly == 2009 | yearListIdx_VisOnly == 2013)
# 
# mergedTrain.set<- mergedSegments[keepDates.train,]
# Train_AcOnly.set <- AcOnlySegments[keepDates_AcOnly.train,]
# Train_VisOnly.set<- VisOnlySegments[keepDatesVisOnly.train,]
# 
# mergedTest.set<- mergedSegments[keepDates.test,]
# Test_AcOnly.set<- AcOnlySegments[keepDates_AcOnly.test,]
# Test_VisOnly.set<- VisOnlySegments[keepDatesVisOnly.test,]
# 

# Cleveland dot plots:
# no transforms
# plot.cleveland(mergedTrain.set,covarList,FALSE,paste0('AcousticAndVisual_',SP))
# plot.cleveland(Train_AcOnly.set,covarList,FALSE,paste0('AcousticOnly_',SP))
# plot.cleveland(Train_VisOnly.set,covarList,FALSE,paste0('VisualOnly_',SP))
# 

# with transformations
# decided to exclude CHL8day (bad distribution), TKE surface current(outliers, plus redundant),
# SST Monthly climate (Looks the same as 8day climate), SSH Monthly climate (same as 8 day climate),
# bathymetry (not normally distributed for fixed sites)
# covarList2 <- c("Density","SST","SSH")

# transformList <- c("none","none","none")
# covarList2 <- c("SST","SSH","CHL","HYCOM_MLD",
#                 "HYCOM_SALIN_0","HYCOM_DIR_0",
#                 "HYCOM_MAG_0",
#                 "HYCOM_UPVEL_50","FrontDist_Cayula",
#                 "EddyDist","Neg_EddyDist","DayOfYear",
#                 "fac1")
# 
# transformList <- c("none","none","log10","log10",
#                    "none","none",
#                    "log10",
#                    "none","log10","none",
#                    "none","none",
#                    "none")
# 
# 
# # restrict covariates again to limited set
# mergedTrain.set2<- mergedTrain.set[,covarList2]
# mergedTest.set2<- mergedTest.set[,covarList2]
# Train_AcOnly.set2<- Train_AcOnly.set[,covarList2]
# Test_AcOnly.set2<- Test_AcOnly.set[,covarList2]
# Train_VisOnly.set2<- Train_VisOnly.set[,covarList2]
# Test_VisOnly.set2<- Test_VisOnly.set[,covarList2]
# 
# 
# transformedCovars.train <- transform.covars(mergedTrain.set2,covarList2,transformList)
# transformedCovars.test <- transform.covars(mergedTest.set2,covarList2,transformList)
# 
# transformedCovars_AcOnly.train <- transform.covars(Train_AcOnly.set2,covarList2,transformList)
# transformedCovars_AcOnly.test <- transform.covars(Test_AcOnly.set2,covarList2,transformList)
# transformedCovars_VisOnly.train <- transform.covars(Train_VisOnly.set2,covarList2,transformList)
# transformedCovars_VisOnly.test <- transform.covars(Test_VisOnly.set2,covarList2,transformList)
# 
# plotCols = colnames(transformedCovars.train)[1:length(covarList2)-1]
# plot.cleveland(transformedCovars.train,plotCols,TRUE,paste0('AcousticAndVisual_',SP))
# plot.cleveland(transformedCovars_AcOnly.train,plotCols,TRUE,paste0('AcousticOnly_',SP))
# plot.cleveland(transformedCovars_VisOnly.train,plotCols,TRUE,paste0('VisualOnly_',SP))
# 
# # presence absence histograms
# plot.covarDensity(transformedCovars.train,plotCols,mergedTrain.set$Density,paste0('AcousticAndVisual_',SP))
# plot.covarDensity(transformedCovars_AcOnly.train,plotCols,Train_AcOnly.set$Density,paste0('AcousticOnly_',SP))
# plot.covarDensity(transformedCovars_VisOnly.train,plotCols,Train_VisOnly.set$Density,paste0('VisualOnly_',SP))
# 
# # correlation
# # without transform
# png(paste(outDir,SP,'_correlations_noTransform.png',sep=''), width = 2000, height = 1600)
# pairs.panels(mergedTrain.set2[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off()
# png(paste(outDir,SP,'_correlations_noTransform_AcOnly.png',sep=''), width = 2000, height = 1600)
# pairs.panels(Train_AcOnly.set2[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off()
# png(paste(outDir,SP,'_correlations_noTransform_visOnly.png',sep=''), width = 2000, height = 1600)
# pairs.panels(Train_VisOnly.set2[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off()
# 
# # with transform
# png(paste(outDir,SP,'_correlations_withTransform.png',sep=''), width = 2000, height = 1600)
# pairs.panels(transformedCovars.train[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off() 
# png(paste(outDir,SP,'_correlations_withTransform_AcOnly.png',sep=''), width = 2000, height = 1600)
# pairs.panels(transformedCovars_AcOnly.train[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off() 
# png(paste(outDir,SP,'_correlations_withTransform_visOnly.png',sep=''), width = 2000, height = 1600)
# pairs.panels(transformedCovars_VisOnly.train[,1:length(covarList2)-1], ellipses=FALSE, method = "spearman",cex.cor=.75)
# dev.off() 
# cat("Exploratory plots done\n")
# 
# 
# save(transformedCovars.train,transformedCovars.test,
#      mergedTest.set,mergedTrain.set,
#      transformedCovars_AcOnly.train,transformedCovars_AcOnly.test,
#      transformedCovars_VisOnly.train,transformedCovars_VisOnly.test,
#      Train_AcOnly.set,Train_VisOnly.set,Test_AcOnly.set,Test_VisOnly.set,
#      file = paste(outDir,SP,'TrainAndTestSets.Rdata',sep=''))

# ############# Other Recipes #########
# 
# # - To keep data needed for map predictions, add the following cotrol item to gamm call
# #       control = list(keepData=TRUE)
# # - After running the gam: 
# # model <- encounterAll_gam
# # coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
# # modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, y, NULL, NULL, coordinateSystem, model)
# # - Then save 'model' to file (.Rdata)
# # save(model, modelMetadata, file = paste(outDir,SP,'_GAM_VisOnly_TF_pruned.Rdata',sep='')) 
# 
# 
# # - To output model summary text to file
# # sink(paste(outDir,SP,'_GAMM_density_full.txt'))
# # summary(encounterAll$gam)
# # gam.check(encounterAll$gam)
# # sink()
# 
# 
# # - To calculate, plot and save residuals to text file
# # rsd <-residuals.gam(encounterAll)
# # png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
# # plot(rsd)
# # dev.off() 
# # 
# # - To plot gam smooths on fewer pages:
# # plot(encounterAll$gam,pages=1,ylim =c(-2,2))
# 
