# Acoustovisual density estimation
library(mrds)
library(psych)
library(lubridate)
library(magic)
library(mgcv)
library(openair)
library(prodlim)

setwd("E:/NASData")
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_missingdata.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_cleveland.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/transform_covars.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_covarDensity.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/GetModelMetadata.R')

#### Parameters needed:  ####
outDir <- "E:/NASData/ModelData/" 
## Detection files
SP <- "Zc" # "Ssp"

acousticSegFile <- "E:/NASData/ALLSITES_segments_daily.csv" # acoustic input file
acousticDensityFile <-"E:/NASData/MC_GC_DT_binsize000800_Group_density_Cuviers.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
# To load this, use: surveyArea <- readShapeSpatial(surveyAreaFile)

## Species/Platform/model info
# Species Category (used for file naming)

# Species names
# Visual Codes
SPC_vis <- c("Cuvier's beaked whale", "unid. Ziphiid","unid. Mesoplodont")# "Gervais' beaked whale", "Beaked Whale","unid. Mesoplodont"
#SPC_vis <- c("Atlantic spotted dolphin", "Striped dolphin","Pantropical spotted dolphin",
#             "Spinner dolphin","Stenella sp.","Clymene dolphin")
# SPC_vis <- c("Risso's dolphin")

# Platform Codes (visual only)
PLC <- c("GU","OR")


# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- FALSE # can be true or false

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# #Vis_TruncDist_allBW.Rdata" # <- With spares data, you could produce a visual detection function using all beaked whales, 
                                        # then  use that here, to only estimate habitat model for Ziphiid, for instance.
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- TRUE  # set to true if you need to match density estimates with environmental parameters
  
visG0 <- .27# mean(c(.42,.37)) #= dolphin #Beaked whale=.27 # from palka 2006 table 5, 2004 survey

############### Begin Calculations ################
graphics.off()
closeAllConnections()

## Load & Prune Acoustic data
cat("Loading acoustic data \n")
if (matchACSegs){
  # Load acoustic segments and densities
  acSegmentsAll <- read.csv(acousticSegFile, header = TRUE,na.strings=c(""," ","NA","-99999","-9999","NaN"))
  acDensityAll <- read.csv(acousticDensityFile, header = TRUE,na.strings=c(""," ","NA","-99999","-9999","NaN"))
  acSegmentsAll$XLSDATE = as.Date(acSegmentsAll$XLSDATE,"%m/%d/%Y")
  acDensityAll$xlsDate = as.Date(acDensityAll$xlsDate,"%m/%d/%Y")#"%m/%d/%Y %H:%M"
  nCol <- length(colnames(acDensityAll))
  
  keepPoints <- which(acDensityAll$xlsDate >= "2011-01-01" & acDensityAll$xlsDate < "2014-01-01")
  acDensityAll <- acDensityAll[keepPoints,]

  # Match segments to density datapoints
  covarNames = names(acSegmentsAll[5:length(names(acSegmentsAll))])
  acDensityAll[,covarNames] <- NA
  for (iR in 1:nrow(acDensityAll)){
    # find all the segments with matching latitudes
    goodLat <- which(acSegmentsAll$LAT %in% acDensityAll$Lat[iR])
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
}else {
  load(acDensityFile)
}
# Exclude partial weeks, and extract the right density estimate type (cue or group)
#fullWeeks <- which(acSegmentsAll$PartialWeek == 0)
#acSegmentsFull <- acSegmentsAll#[fullWeeks,]
rm(acSegmentsAll)
  
###########################

## Load Visual data
cat("Loading visual data\n")

load(visDataFile)
visSegments <- read.csv(visSegmentsFile, header = TRUE,na.strings=c(""," ","NA","-99999.0000","-99999"))
visSegments$date <- as.Date(visSegments$date,"%m/%d/%Y")
## Process Visual data to determine detection probabilities and strip widths
visSpIdx <- NULL
for (spname in SPC_vis){
  visSpIdx <- c(visSpIdx,which(grepl(spname, visData$commonname)))
}

# prune out off-effort sightings
spIdxON <- visSpIdx[which(as.logical(visData$effort[visSpIdx]))]

# prune out sightings with an angle over 90 deg
spIdxON2 <- spIdxON[which(visData$relbear[spIdxON]<90)]

# Populate truncated column of on effort sightings of species of interest with zeros
visData$Truncated <- 1
visData$Truncated[spIdxON2] <- 0


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
    PLCspIdxOn <- spIdxON2[which(grepl(i,visData$ship[spIdxON]))] 
    
    # Make dataframe with the inputs that the ddf distance function wants
    ddfData$observer <- rep(1,length(PLCspIdxOn))
    ddfData$detected <- rep(1,length(PLCspIdxOn))
    ddfData$object <- (1:length(PLCspIdxOn))
    ddfData$distance <- visData$transect_distm[PLCspIdxOn]
    ddfData$size <- visData$size[PLCspIdxOn]
    ddfData$seastate <- visData$seastate[PLCspIdxOn]
    ddfData$swell <- visData$swell[PLCspIdxOn]
    ddfData$vis <- visData$vis[PLCspIdxOn]
    
    
    # Compute untruncated detection function
    cat("Calculating basic fit with non-truncated data and half-normal key, no covariates.\n")
    
    detFun_noTrunc <-ddf(method = 'ds', dsmodel =~ mcds(key = 'hn', formula = ~ 1),
                         data = as.data.frame(ddfData), meta.data = list(binned=F,left=0),
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
                  data = as.data.frame(ddfData), meta.data = list(binned=F, width=tDist[nPlatform], left=0))
#         df <- ds(as.data.frame(ddfData), truncation = tDist[nPlatform], order = NULL, transect = "line", key = keyListInit[[i1]],
#                   monotonicity = "weak")
        
      } else {
        df <-ddf(method ='ds', dsmodel =~ mcds(key = keyListInit[[i1]], formula = ~ 1,
                adj.series = adjInit[[i1]], adj.order = adjOrderInit[[i1]]), data = as.data.frame(ddfData),
                meta.data = list(binned=F, width=tDist[nPlatform],left=0))
#         df <- ds(as.data.frame(ddfData), truncation = tDist[nPlatform], transect = "line", key = keyListInit[[i1]],
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

    cat("Done fitting models")
    
    # Put all combinations together, and see which one has the lowest AIC
    aicList<-c(aicList1)#,aicList2) (commented out part associated with covariate-models
    keyList <- c(keyList1)#,keyList2)
    adjList <- c(adjStr)#,cSetStr)
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
##

# Put effective strip widths calculated above into Segments table
visSegments["ESW"] <- 0
for (iP in 1:length(PLC)){
  thisSegList <-which(visSegments$ship==PLC[iP])
  visSegments$ESW[thisSegList] <- predict(detFunByPlatform[[iP]],esw=TRUE)$fitted[1]
}

# get rid of off effort segments
visSeg_OnEffort <- visSegments[visSegments$effort==1,]

# Tally encounters by segment
prunedSightings <- visData[visData$Truncated==0,] # get all of the non-truncated sightings

# match assign segment to each sighting
for (iSight in 1:length(prunedSightings$date)){
  sightDate <- prunedSightings$date[iSight]
  onThisDay <- which(visSeg_OnEffort$date == sightDate)
  if (length(onThisDay)>0) {
    minIdx <- which.min(abs(visSeg_OnEffort$transect[onThisDay]-prunedSightings$transect[iSight]))
    prunedSightings$Segment[iSight] <-onThisDay[minIdx]
  }else { # handle case where there is no match (why would this happen?)
    cat(paste("Warning: Missing effort segment for sighting on", sightDate,"\n"))
    prunedSightings$Segment[iSight] <- NaN
  }
}

segTally <- as.data.frame(table(prunedSightings$Segment)) # this gives you a list of segments containing sightings

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
for (iSeg in segTally[,1]){
  
  # determine the row number of all sighting rows matching this segment number
  segIdx <- which(prunedSightings$Segment == iSeg)
  SegObjID_Idx <- which(visSeg_OnEffort$OBJECTID == iSeg)
  
  if (length(SegObjID_Idx)>0){
    visSeg_OnEffort$sp_count[SegObjID_Idx] <- sum(prunedSightings$size[segIdx])
    if (visSeg_OnEffort$sp_count[SegObjID_Idx] >0) {# this should always be true...
      visSeg_OnEffort$sp_present[SegObjID_Idx] <- 1
    }
  }
}

# account for G0 in encounters
visSeg_OnEffort$sp_count_g0adj <- visSeg_OnEffort$sp_count/visG0

# Estimate surveyed area
visSeg_OnEffort$EffectiveArea <- (2*visSeg_OnEffort$ESW/1000)*(visSeg_OnEffort$SegmentLength/1000)
visSeg_OnEffort$Density <- visSeg_OnEffort$sp_count_g0adj/visSeg_OnEffort$EffectiveArea


###################################
## Visual and Acoustic

AcOnlySegments <- NULL
yearListTemp <- as.numeric(format(acDensityAll$xlsDate,"%Y"))
siteYear <- NULL
nAc <- length(acDensityAll$HYCOM_MAG_100)
siteYear$Year <-yearListTemp
siteYear$Site <-acDensityAll$Site
siteYear <- as.data.frame(siteYear)
uSiteYear <- unique((siteYear))
AcOnlySegments$Density <- rep(NA,times = nAc)
for (uR in 1:nrow(uSiteYear)){
  # Normalize acoustic density estimated by deployment
  thisSet <- which(as.logical(row.match(siteYear,uSiteYear[uR,])))
  thisSet_gt0 <- which(acDensityAll$meanDensity[thisSet]>0)
  quant95 <-quantile(acDensityAll$meanDensity[thisSet[thisSet_gt0]],probs = .95,na.rm = TRUE)
  AcOnlySegments$Density[thisSet] <- acDensityAll$meanDensity[thisSet]/quant95
}

AcOnlySegments$date <- as.Date(acDensityAll$xlsDate,"%m/%d/%Y") #date
AcOnlySegments$lat <- acDensityAll$Lat
AcOnlySegments$long <- acDensityAll$Long
AcOnlySegments$Category<- rep(2,length(acDensityAll$Long))
AcOnlySegments$SST <- acDensityAll$SST_DAILY_CMC.L4
AcOnlySegments$SSH <- acDensityAll$SSH_DAILY_AVISO
AcOnlySegments$CHL <- acDensityAll$CHL_8DAY_NASA
AcOnlySegments$HYCOM_QTOT <- acDensityAll$HYCOM_QTOT
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
AcOnlySegments$EddyDist <- acDensityAll$EDDYDIST
AcOnlySegments$Type <- rep(2,times = nAc)

VisOnlySegments <- NULL
VisOnlySegments$date <- as.Date(visSeg_OnEffort$date_Converted,"%m/%d/%Y") #date
VisOnlySegments$lat <- visSeg_OnEffort$Lat
VisOnlySegments$long <- visSeg_OnEffort$Long
VisOnlySegments$Category<- rep(1,length(visSeg_OnEffort$Long))
# mergedSegments$ESW <- c(acDensityAll$BW_ESW)
VisOnlySegments$Density <- visSeg_OnEffort$Density/
  quantile(visSeg_OnEffort$Density[which(visSeg_OnEffort$Density>0)],probs = .95,na.rm = TRUE)
VisOnlySegments$SST <- visSeg_OnEffort$SST_daily_CMC_L4_GLOB
VisOnlySegments$SSH <- visSeg_OnEffort$SSH_daily_aviso_double
VisOnlySegments$CHL <- visSeg_OnEffort$CHl_8Day_NASA
VisOnlySegments$HYCOM_QTOT <- visSeg_OnEffort$HYCOM_qTot
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
nVis <- length(visSeg_OnEffort$HYCOM_upVel_100)
VisOnlySegments$Type <- rep(1,times = nVis)

cat("Merging Visual and Acoustic Segments\n")

# Merge visual and acoustic segments into one big dataframe
mergedSegments <- NULL
mergedSegments$date <- c(as.Date(visSeg_OnEffort$date_Converted,"%m/%d/%Y"),as.Date(acDensityAll$xlsDate,"%m/%d/%Y")) #date
mergedSegments$lat <- c(visSeg_OnEffort$Lat,acDensityAll$Lat)
mergedSegments$long <- c(visSeg_OnEffort$Long,acDensityAll$Long)
mergedSegments$Category<- c(rep(1,length(visSeg_OnEffort$Long)),rep(2,length(acDensityAll$Long)))
# mergedSegments$ESW <- c(acDensityAll$BW_ESW)
mergedSegments$Density <- c(VisOnlySegments$Density,
                            AcOnlySegments$Density)
mergedSegments$SST <- c(visSeg_OnEffort$SST_daily_CMC_L4_GLOB,acDensityAll$SST_DAILY_CMC.L4)
mergedSegments$SSH <- c(visSeg_OnEffort$SSH_daily_aviso_double,acDensityAll$SSH_DAILY_AVISO)
mergedSegments$CHL <- c(visSeg_OnEffort$CHl_8Day_NASA, acDensityAll$CHL_8DAY_NASA)
mergedSegments$HYCOM_QTOT <- c(visSeg_OnEffort$HYCOM_qTot, acDensityAll$HYCOM_QTOT)
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
mergedSegments$EddyDist <- c(visSeg_OnEffort$EddyDist,acDensityAll$EDDYDIST)


mergedSegments$Type <- c(rep(1,times = nVis),rep(2,times = nAc))


mergedSegments <- as.data.frame(mergedSegments)
AcOnlySegments <- as.data.frame(AcOnlySegments)
VisOnlySegments <- as.data.frame(VisOnlySegments)

covarList<-names(mergedSegments[c(2,5:length(mergedSegments))])
#covarList<- c("Bathymetry","SST_daily","CHl_8Day","HYCOM_dir_daily","HYCOM_mag_daily","HYCOM_wVel_daily",
#              "SSH_daily","CHL_daily","TKE_surfaceCurrent_5day","HYCOM_mld_daily")
###################################
# Oceanographic variables

# Explore the data, graphical output
cat("Begin exploratory plot generation\n")  
# histograms of missing data
percFilled <- plot.missingdata(mergedSegments,covarList,paste0('AcousticAndVisual_',SP)) 
percFilled <- plot.missingdata(AcOnlySegments,covarList,paste0('AcousticOnly_',SP)) 
percFilled <- plot.missingdata(VisOnlySegments,covarList,paste0('VisualOnly_',SP)) 

# Identify and clear problematic outliers
outlierList <-which(mergedSegments$CHL<  -10)
mergedSegments$CHL[outlierList] <- NaN 
outlierList <-which(mergedSegments$FrontDist_Cayula>800000)
mergedSegments$FrontDist_Cayula[outlierList] <- NaN 
outlierList <-which(mergedSegments$Density>10)
mergedSegments$Density[outlierList] <- NaN 

outlierList <-which(AcOnlySegments$CHL<  -10)
AcOnlySegments$CHL[outlierList] <- NaN 
outlierList <-which(AcOnlySegments$FrontDist_Cayula>800000)
AcOnlySegments$FrontDist_Cayula[outlierList] <- NaN 
outlierList <-which(AcOnlySegments$Density>10)
AcOnlySegments$Density[outlierList] <- NaN 

outlierList <-which(VisOnlySegments$CHL<  -10)
VisOnlySegments$CHL[outlierList] <- NaN 
outlierList <-which(VisOnlySegments$FrontDist_Cayula>800000)
VisOnlySegments$FrontDist_Cayula[outlierList] <- NaN 
outlierList <-which(VisOnlySegments$Density>10)
VisOnlySegments$Density[outlierList] <- NaN 

# If you decide from the missing data plots that you want to restrict years going forward:
yearListIdx = as.numeric(format(mergedSegments$date,"%Y"))
yearListIdx_AcOnly = as.numeric(format(AcOnlySegments$date,"%Y"))
yearListIdx_VisOnly = as.numeric(format(VisOnlySegments$date,"%Y"))

keepDates.train <- which(yearListIdx != 2009 & yearListIdx >= 2003 & yearListIdx <= 2012)
keepDates.test <- which(yearListIdx == 2009 | yearListIdx == 2013)
keepDates_AcOnly.train <- which(yearListIdx_AcOnly != 2009 & yearListIdx_AcOnly >= 2003 & yearListIdx_AcOnly <= 2012)
keepDates_AcOnly.test <- which(yearListIdx_AcOnly == 2009 | yearListIdx_AcOnly == 2013)
keepDatesVisOnly.train <- which(yearListIdx_VisOnly != 2009 & yearListIdx_VisOnly >= 2003 & yearListIdx_VisOnly <= 2012)
keepDatesVisOnly.test <- which(yearListIdx_VisOnly == 2009 | yearListIdx_VisOnly == 2013)




mergedTrain.set<- mergedSegments[keepDates.train,]
Train_AcOnly.set<- AcOnlySegments[keepDates_AcOnly.train,]
Train_VisOnly.set<- VisOnlySegments[keepDatesVisOnly.train,]

mergedTest.set<- mergedSegments[keepDates.test,]
Test_AcOnly.set<- AcOnlySegments[keepDates_AcOnly.test,]
Test_VisOnly.set<- VisOnlySegments[keepDatesVisOnly.test,]


# Cleveland dot plots:
# no transforms
plot.cleveland(mergedTrain.set,covarList,FALSE,paste0('AcousticAndVisual_',SP))
plot.cleveland(Train_AcOnly.set,covarList,FALSE,paste0('AcousticOnly_',SP))
plot.cleveland(Train_VisOnly.set,covarList,FALSE,paste0('VisualOnly_',SP))


# with transformations
# decided to exclude CHL8day (bad distribution), TKE surface current(outliers, plus redundant),
# SST Monthly climate (Looks the same as 8day climate), SSH Monthly climate (same as 8 day climate),
# bathymetry (not normally distributed for fixed sites)
# covarList2 <- c("Density","SST","SSH")

# transformList <- c("none","none","none")
covarList2 <- c("SST","SSH","CHL","HYCOM_MLD",
                "HYCOM_SALIN_0",
                "HYCOM_MAG_0",
                "HYCOM_UPVEL_50","FrontDist_Cayula","EddyDist")

transformList <- c("none","none","log10","log10",
                   "none",
                   "log10",
                   "none","log10","none")


# restrict covariates again to limited set
mergedTrain.set2<- mergedTrain.set[,covarList2]
mergedTest.set2<- mergedTest.set[,covarList2]
Train_AcOnly.set2<- Train_AcOnly.set[,covarList2]
Test_AcOnly.set2<- Test_AcOnly.set[,covarList2]
Train_VisOnly.set2<- Train_VisOnly.set[,covarList2]
Test_VisOnly.set2<- Test_VisOnly.set[,covarList2]



transformedCovars.train <- transform.covars(mergedTrain.set2,covarList2,transformList)
transformedCovars.test <- transform.covars(mergedTest.set2,covarList2,transformList)

transformedCovars_AcOnly.train <- transform.covars(Train_AcOnly.set2,covarList2,transformList)
transformedCovars_AcOnly.test <- transform.covars(Test_AcOnly.set2,covarList2,transformList)
transformedCovars_VisOnly.train <- transform.covars(Train_VisOnly.set2,covarList2,transformList)
transformedCovars_VisOnly.test <- transform.covars(Test_VisOnly.set2,covarList2,transformList)

plot.cleveland(transformedCovars.train,colnames(transformedCovars.train),TRUE,paste0('AcousticAndVisual_',SP))
plot.cleveland(transformedCovars_AcOnly.train,colnames(transformedCovars.train),TRUE,paste0('AcousticOnly_',SP))
plot.cleveland(transformedCovars_VisOnly.train,colnames(transformedCovars.train),TRUE,paste0('VisualOnly_',SP))

# presence absence histograms
plot.covarDensity(transformedCovars.train,colnames(transformedCovars.train),mergedTrain.set$Density,paste0('AcousticAndVisual_',SP))
plot.covarDensity(transformedCovars_AcOnly.train,colnames(transformedCovars_AcOnly.train),Train_AcOnly.set$Density,paste0('AcousticOnly_',SP))
plot.covarDensity(transformedCovars_VisOnly.train,colnames(transformedCovars_VisOnly.train),Train_VisOnly.set$Density,paste0('VisualOnly_',SP))

# correlation
# without transform
png(paste(outDir,SP,'_correlations_noTransform.png',sep=''), width = 2000, height = 1600)
pairs.panels(mergedTrain.set2, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off()
png(paste(outDir,SP,'_correlations_noTransform_AcOnly.png',sep=''), width = 2000, height = 1600)
pairs.panels(Train_AcOnly.set2, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off()
png(paste(outDir,SP,'_correlations_noTransform_visOnly.png',sep=''), width = 2000, height = 1600)
pairs.panels(Train_VisOnly.set2, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off()

# with transform
png(paste(outDir,SP,'_correlations_withTransform.png',sep=''), width = 2000, height = 1600)
pairs.panels(transformedCovars.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off() 
png(paste(outDir,SP,'_correlations_withTransform_AcOnly.png',sep=''), width = 2000, height = 1600)
pairs.panels(transformedCovars_AcOnly.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off() 
png(paste(outDir,SP,'_correlations_withTransform_visOnly.png',sep=''), width = 2000, height = 1600)
pairs.panels(transformedCovars_VisOnly.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off() 
cat("Exploratory plots done\n")

###########################
# Run & evaluate models

# # Presence absence
# yBinomial <- mergedTrain.set$SpPresent
# 
# cat("Run full GAM on presence absence data with shrinkage\n")

yAcOnly_TF <- as.logical(Train_AcOnly.set$Density>0)
# myts_AcOnly = ts(Train_AcOnly.set$Density[which(Train_AcOnly.set$lat == 28.84625)],start = 1, frequency = 1)
# tsdiag(arima(myts_AcOnly))

cat("Run full binomial GAM on Acoustic only data with shrinkage\n")#random = list(fac1=~1),
gam_full_AcOnly_TF <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                           + s(SSH, bs="ts", k=kVal)
                           + s(log10_FrontDist_Cayula, bs="ts", k=kVal)
                           + s(Neg_EddyDist, bs="ts", k=kVal)
                           + s(DayOfYear, bs="cp", k=kVal)
                           + s(log10_HYCOM_MAG_100,bs="ts", k=kVal)
                           + s(HYCOM_SALIN_100, bs="ts", k=kVal),
                           data = transformedCovars_AcOnly.train,
                           na.action = na.omit,family = quasibinomial(),
                           correlation = corAR1(form=~1|fac1))#


# Save summary to text file
sink(paste(outDir,SP,'_GAMM_presence_full.txt'))
summary(gam_full_AcOnly_TF$gam)
sink()
# 
# # Calculate and save residuals to text file
# rsd <-residuals.gam(presAbsGAMAll)
# png(paste(outDir,SP,'_residuals_presence_full.png',sep=''), width = 1000, height = 800)
# plot(rsd)
# dev.off() 

# # Density
# kVal = 8
# yAcOnly <- (Train_AcOnly.set$Density)
# # myts_AcOnly = ts(Train_AcOnly.set$Density[which(Train_AcOnly.set$lat == 28.84625)],start = 1, frequency = 1)
# # tsdiag(arima(myts_AcOnly))
# cat("Run full GAM on Acoustic only data with shrinkage\n")#correlation = corAR1(form=~1|mergedTrain.set$Category)
# gam_full_AcOnly <- gam(yAcOnly~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
#                     + s(log10_HYCOM_MLD, bs="ts",k=kVal)+ s(HYCOM_SALIN_0, bs="ts",k=kVal)
#                     + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)+ s(HYCOM_UPVEL_50, bs="ts",k=kVal)
#                     + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
#                     data = transformedCovars_AcOnly.train,  method = "GCV.Cp",
#                     na.action = na.omit,family=tw())#
# 
# sink(paste(outDir,SP,'_GAM_full_AcOnly.txt'))
# summary(gam_full_AcOnly)
# gam.check(gam_full_AcOnly)
# sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(gam_full_AcOnly)
png(paste(outDir,SP,'_residuals_density_AcOnly_Full.png', sep=''), width = 1000, height = 800)
plot(rsd)
dev.off()


yAcOnly_TF <- as.logical(Train_AcOnly.set$Density>0)
# myts_AcOnly = ts(Train_AcOnly.set$Density[which(Train_AcOnly.set$lat == 28.84625)],start = 1, frequency = 1)
# tsdiag(arima(myts_AcOnly))
cat("Run full binomial GAM on Acoustic only data with shrinkage\n")#correlation = corAR1(form=~1|mergedTrain.set$Category)
gam_full_AcOnly_TF <- gam(yAcOnly_TF~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                       + s(HYCOM_SALIN_0, bs="ts",k=kVal)
                       + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)
                       + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
                       data = transformedCovars_AcOnly.train,  method = "GCV.Cp",
                       na.action = na.omit,family=binomial(),control = list(keepData=TRUE))#
model <- gam_full_AcOnly_TF
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars_AcOnly.train, NULL, yAcOnly_TF, NULL, NULL, coordinateSystem, model)
save(model, modelMetadata, file = paste(outDir,SP,'_GAM_AcOnly_TF_pruned.Rdata',sep='')) 
sink(paste(outDir,SP,'_GAM_full_AcOnly_binomial.txt'))
summary(gam_full_AcOnly_TF)
gam.check(gam_full_AcOnly_TF)
sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(gam_full_AcOnly_TF)
png(paste(outDir,SP,'_residuals_binom_AcOnly_Full.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 


yVisOnly <- (Train_VisOnly.set$Density)
plot(yVisOnly)
cat("Run full GAM on Visual only data with shrinkage\n")#correlation = corAR1(form=~1|mergedTrain.set$Category)
gam_full_VisOnly <- gam(yVisOnly~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                       + s(log10_HYCOM_MLD, bs="ts",k=kVal)+ s(HYCOM_SALIN_0, bs="ts",k=kVal)
                       + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)+ s(HYCOM_UPVEL_50, bs="ts",k=kVal)
                       + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
                       data = transformedCovars_VisOnly.train,  method = "GCV.Cp",
                       na.action = na.omit,family=tw())#

sink(paste(outDir,SP,'_GAM_full_VisOnly.txt'))
summary(gam_full_VisOnly)
gam.check(gam_full_VisOnly)
sink()

yVisOnly_TF <- (Train_VisOnly.set$Density>0)
plot(yVisOnly_TF)
kVal<-8
cat("Run full binomial GAM on Visual only data with shrinkage\n")#correlation = corAR1(form=~1|mergedTrain.set$Category)
gam_full_VisOnly_TF <- gam(yVisOnly_TF~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts",k=kVal)
                        + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)
                        + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
                        data = transformedCovars_VisOnly.train,  method = "GCV.Cp",
                        na.action = na.omit,family=binomial(),control = list(keepData=TRUE))#
model <- gam_full_VisOnly_TF
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars_VisOnly.train, NULL, yVisOnly_TF, NULL, NULL, coordinateSystem, model)
save(model, modelMetadata, file = paste(outDir,SP,'_GAM_VisOnly_TF_pruned.Rdata',sep='')) 


sink(paste(outDir,SP,'_GAM_full_VisOnly_binomial.txt'))
summary(gam_full_VisOnly_TF)
gam.check(gam_full_VisOnly_TF)
sink()
plot(gam_full_VisOnly_TF,ylim = c(-5,5),pages = 1)


yTF <- (mergedTrain.set$Density>0)
plot(yTF)
cat("Run full TF GAMM on density data with shrinkage\n")#correlation = corAR1(form=~1|mergedTrain.set$Category)
myCat <- as.factor(mergedTrain.set$Category)
myWeights <- rep(1,times = length(myCat))
AcTrainSize <- length(which(myCat==2))
VisTrainSize <- length(which(myCat==1))
myWeights[which(myCat==2)] <- VisTrainSize/AcTrainSize

Sub = transformedCovars.train
encounterAll_gamm_TF <- gamm(yTF ~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                     + s(HYCOM_SALIN_0, bs="ts",k=kVal)
                     + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)
                     + s(EddyDist, bs="ts",k=kVal),
                     data = transformedCovars.train, random=list(myCat=~1),weights = myWeights, method = "GCV.Cp",
                     na.action = na.omit,family=binomial(),niterPQL = 40,control = list(keepData=TRUE,maxIter=200,niterEM=0))#,family= Tweedie(p=1.4)

sink(paste(outDir,SP,'_GAMM_full_TF.txt'))
summary(encounterAll_gamm_TF$gam)
gam.check(encounterAll_gamm_TF$gam)
sink()
# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll_gamm_TF$gam)
model <- encounterAll_gamm_TF$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, yTF, NULL, NULL, coordinateSystem, model)
save(model, modelMetadata, file = paste(outDir,SP,'_GAMM_TF_pruned.Rdata',sep='')) 


y <- (mergedTrain.set$Density)
encounterAll_gamm <- gamm(y~ s(SST, bs="ts",k=5)+ s(SSH, bs="ts",k=5)+ s(log10_CHL, bs="ts",k=5)
                             + s(log10_HYCOM_MLD, bs="ts",k=5)+ + s(HYCOM_DIR_0, bs="cc",k=5) + s(HYCOM_SALIN_0, bs="ts",k=5)
                             + s(log10_HYCOM_MAG_0, bs="ts",k=5)+ s(HYCOM_UPVEL_100, bs="ts",k=5),
                             data = transformedCovars.train, random=list(myCat=~1),weights = myWeights,method = "GCV.Cp",
                             na.action = na.omit,family= Tweedie(p=1.2))#

sink(paste(outDir,SP,'_GAMM_full.txt'))
summary(encounterAll_gamm$gam)
gam.check(encounterAll_gamm$gam)
sink()
# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll_gamm$gam)


encounterAll_gam_TF <- gam(yTF~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                     + s(log10_HYCOM_MLD, bs="ts",k=kVal)+ s(HYCOM_SALIN_0, bs="ts",k=kVal)
                     + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)+ s(HYCOM_UPVEL_50, bs="ts",k=kVal)
                     + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
                     data =transformedCovars.train,  method = "GCV.Cp",na.action = na.omit,family=binomial(), weights = myWeights,
                     control = list(keepData=TRUE))#

sink(paste(outDir,SP,'_GAM_full_TF.txt'))
summary(encounterAll_gam_TF)
gam.check(encounterAll_gam_TF)
sink()


encounterAll_gam <- gam(y~ s(SST, bs="ts",k=kVal)+ s(SSH, bs="ts",k=kVal)+ s(log10_CHL, bs="ts",k=kVal)
                         + s(log10_HYCOM_MLD, bs="ts",k=kVal)+ s(HYCOM_SALIN_0, bs="ts",k=kVal)
                         + s(log10_HYCOM_MAG_0, bs="ts",k=kVal)+ s(HYCOM_UPVEL_50, bs="ts",k=kVal)
                         + s(log10_FrontDist_Cayula, bs="ts",k=kVal)+ s(EddyDist, bs="ts",k=kVal),
                         data =Sub,  method = "GCV.Cp",na.action = na.omit,family=tw(), weights = myWeights,
                         control = list(keepData=TRUE))#

sink(paste(outDir,SP,'_GAM_full_TF.txt'))
summary(encounterAll_gam)
gam.check(encounterAll_gam)
sink()
model <- encounterAll_gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, y, NULL, NULL, coordinateSystem, model)
# save model
save(model, modelMetadata, file = paste(outDir,SP,'_GAM_pruned.Rdata',sep=''))


# png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
# plot(rsd)
# dev.off() 


# encounterAll <- gam(y ~ s(SST, bs="ts",k=5) + s(SSH, bs="ts",k=5) + s(log10_CHL, bs="ts",k=5) + s(log10_HYCOM_MLD, bs="ts",k=5) + s(HYCOM_EMP, bs="ts",k=5) 
#          + s(HYCOM_DIR_0, bs="cc",k=5) + s(HYCOM_DIR_100, bs="cc",k=5)
#          + s(HYCOM_SALIN_0, bs="ts",k=5) + s(HYCOM_SALIN_100, bs="ts",k=5) + s(HYCOM_SALIN_800, bs="ts",k=5)
#          + s(log10_HYCOM_MAG_0, bs="ts",k=5) + s(log10_HYCOM_MAG_100, bs="ts",k=5) + s(HYCOM_MAG_800, bs="ts",k=5)
#          + s(HYCOM_UPVEL_100, bs="ts",k=5) + s(HYCOM_UPVEL_800),
#          method = "REML", data = transformedCovars.train, family = tw(),offset = log(mergedTrain.set$EffectiveArea),
#          na.action = na.omit)# 

# Output summary text to file
sink(paste(outDir,SP,'_GAMM_density_full.txt'))
summary(encounterAll$gam)
gam.check(encounterAll$gam)
sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll)
png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 

plot(encounterAll$gam,pages=1,ylim =c(-2,2))

cat("Run reduced GAMM without unused variables\n")
# look at that output, some covariates have been shrunk down to nothing, so remove them
ctl <- gam.control()
ctl$keepData = TRUE
# encounterAll <- gamm(y ~ s(SSH, bs="ts",k=5)+ s(log10_HYCOM_MLD, bs="ts", k=5)
#                     + s(HYCOM_DIR_0, bs="cc", k=5)
#                     + s(HYCOM_SALIN_0, bs="ts", k=5),
#                     control = list(keepData=TRUE), random=list(myCat=~1),weights = myWeights,
#                     method = "REML", data = Sub,  family = Tweedie(1.4),
#                     na.action = na.omit)#
encounterAll <- gamm(y ~  s(SST, bs="ts",k=kVal) + s(SSH, bs="ts", k=kVal) +
                     +  s(log10_CHL, bs="ts",k=kVal)+ s(log10_HYCOM_MLD, bs="ts",k=kVal) +s(HYCOM_SALIN_0, bs="ts", k=kVal)
                     + s(log10_HYCOM_MAG_0, bs="ts", k=kVal)+ s(HYCOM_UPVEL_50, bs="ts", k=kVal)
                     + s(log10_FrontDist_Cayula, bs="ts", k=kVal)+ s(EddyDist, bs="ts", k=kVal),
                     control = list(keepData = TRUE), random=list(myCat=~1),weights = myWeights,
                     method = "GCV.Cp", data = Sub,  family = Tweedie(1.4),
                     na.action = na.omit)#
 
plot(encounterAll$gam,pages=1,ylim =c(-2,2))

# Output summary text to file
sink(paste(outDir,SP,'_GAM_density_pruned.txt'))
summary(encounterAll$gam)
gam.check(encounterAll$gam)
sink()

plot(encounterAll$gam)

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll$gam)
png(paste(outDir,SP,'_residuals_density_pruned.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 

# save model


# Predict on test data
yTest1 <- mergedTest.set$Density[which(!is.na(mergedTest.set$Density))]
yTest <- mergedTest.set$Density>0
pred <- predict.gam(encounterAll_gamm_TF$gam,transformedCovars.test, type = 'response',na.action = na.omit)
plot(pred)

# Density 
