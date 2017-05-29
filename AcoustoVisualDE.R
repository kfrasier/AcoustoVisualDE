# Acoustovisual density estimation
library(mrds)
library(psych)
library(lubridate)
library(magic)
library(mgcv)
setwd("E:\\NASData")
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_missingdata.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_cleveland.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/transform_covars.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_covarDensity.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/GetModelMetadata.R')

#### Parameters needed:  ####
outDir <- "E:/NASData/ModelData/" 
## Detection files
acousticSegFile <- "E:/NASData/Acoustic_Segments_Ssp_only.csv" # acoustic input file

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
# To load this, use: surveyArea <- readShapeSpatial(surveyAreaFile)

## Species/Platform/model info
# Species Category (used for file naming)
SP <- "Ssp" # "Zc"

# Species names
# Visual Codes
# SPC_vis <- c("Cuvier's beaked whale", "unid. Ziphiid")# "Gervais' beaked whale", "Beaked Whale","unid. Mesoplodont"
SPC_vis <- c("Atlantic spotted dolphin", "Striped dolphin","Pantropical spotted dolphin",
             "Spinner dolphin","Stenella sp.","Clymene dolphin")


# Acoustic Codes
# SPC_ac <- c("Cuvier's beaked whale")
SPC_ac <- c("Ssp")

# Platform Codes (visual only)
PLC <- c("GU","OR")

# Use cue or group-based density estimates?
acModel <- "cue"  # can be 'cue' or 'group'

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- TRUE # can be true or false

# The name of the visual detection function file. 
detFunFile <- "Vis_TruncDist_Ssp_only.Rdata" #Vis_TruncDist_allBW.Rdata" # <- With spares data, you could produce a visual detection function using all beaked whales, 
                                        # then  use that here, to only estimate habitat model for Ziphiid, for instance.
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

visG0 <- mean(c(.42,.37))#Beaked whale:.27 # from palka 2006 table 5, 2004 survey




############### Begin Calculations ################
graphics.off()
closeAllConnections()

## Load & Prune Acoustic data
cat("Loading acoustic data \n")

# Note: I don't think there need to make both "Sightings" and "Segments" equivalents for acoustic
# data, because there is a 1-t0-1 relationship between the two in this case. There is a data point
# for each week, regardless of whether beaked whales were detected.
acSegmentsAll <- read.csv(acousticSegFile, header = TRUE,na.strings=c(""," ","NA","-99999","NaN"))


# Exclude partial weeks, and extract the right density estimate type (cue or group)
fullWeeks <- which(acSegmentsAll$PartialWeek == 0)
acSegmentsFull <- acSegmentsAll[fullWeeks,]

if (acModel == "cue"){
  deIDx <- which(acSegmentsFull$click_de == 1)
} else if (acModel == "group") {
  deIdx <- which(acSegmentsFull$group_de == 1)
} else {
  stop("bad acoustic density estimate type. Options are 'cue' or 'group'.")
}
acSegmentsPruned <- acSegmentsFull[deIDx,]

# Find all of the rows that match the species of interest (can be multiple species names)
acSpIdx <- NULL;
for (spname in SPC_ac){
  acSpIdx <- c(acSpIdx,which(grepl(spname, acSegmentsPruned$commonname)))
}
acSegments <- acSegmentsPruned[acSpIdx,] # prune to retain only those rows


###########################

## Load Visual data
cat("Loading visual data\n")

load(visDataFile)
visSegments <- read.csv(visSegmentsFile, header = TRUE,na.strings=c(""," ","NA","-99999"))
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

#     #### Fitting detection functions with covariates: We can't do this unless we know 
#     # the values of each covariate at all transect segments in order to calculate "Effective Strip Width".
#     # Don't have those values, so we are not doing this.

#     # Iterate over detection functions with covariates and identify AIC for each
#     covarList <- c("size", "seastate") #, 'vis', 'swell')
#     
#     # list of keys: hn, hr
#     keyCovar= c('hn', 'hr')
#     CI <- 1
#     aicList2 <- NULL # store AIC scores 
#     keyList2 <- NULL # store keys scores 
#     cSetStr <- NULL # store the covariate formulas
#     detFun2<- NULL
#     
#     # iterate over the key options
#     for (iKey in keyCovar){
#       
#       # iterate over covariate combinations - using "combn" to come up with different the combinations of covariates
#       for (i2 in 1 : length(covarList)){
#         covarSet <- combn(covarList,i2)
#         nSets <- dim(covarSet)
#         
#         # for each set of covariates, fit a model
#         for (i3 in 1:nSets[2]){
#           cSet <- covarSet[,i3]
#           # nPlus <- nSets[1]-1
#           cSetStr[CI] <- paste('~', paste0(cSet, collapse = " + "))
#           
#           # sometimes models do not converge, use try() to avoid crash if a model fails
#           dF <- NULL
#           try(dF <- ddf(method='ds',dsmodel=~mcds(key=iKey,  formula = cSetStr[CI]),
#                         data = as.data.frame(ddfData),
#                         meta.data = list(binned=F, width=tDist[nPlatform],left=0)))
#           
#           if (is.null(dF)){
#             cat(paste0("Model did not converge: Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
#             
#           }else {
#             detFun2[[CI]] <-dF
#             
#             aicList2[CI] <- detFun2[[CI]]$criterion
#             keyList2[CI] <- iKey
#             cat(paste0("Model result ", CI+i1, ": Key = ",iKey, "; covariates = ",cSetStr[CI],"\n", collapse = ""))
#             cat(paste("AIC =",  round(aicList2[CI], digits=2),"\n", collapse = ""))
#             CI <- CI+1
#           }
#           
#           
#           
#         }
#      }
#     }
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


###################################
## Visual and Acoustic

cat("Merging Visual and Acoustic Segments\n")

# Merge visual and acoustic segments into one big dataframe
mergedSegments <- NULL
mergedSegments$date <- c(as.Date(acSegmentsPruned$date_Converted,"%m/%d/%Y")) #date
mergedSegments$lat <- c(acSegmentsPruned$Lat)
mergedSegments$long <- c(acSegmentsPruned$Long)
mergedSegments$ESW <- c(acSegmentsPruned$BW_ESW)
mergedSegments$SpPresent <- c(acSegmentsPruned$BW_Present)
mergedSegments$SpEncounter <- c(acSegmentsPruned$BW_Encounter)
mergedSegments$SpEncounter_g0adj <- c(acSegmentsPruned$BW_Encounter)
mergedSegments$EffectiveArea <- c(acSegmentsPruned$BW_ESW)


mergedSegments$Bathymetry <- c(acSegmentsPruned$Bathymetry)
mergedSegments$SST_daily <- c(acSegmentsPruned$SST_daily)
mergedSegments$SSH_daily <- c(acSegmentsPruned$SSH_daily)
mergedSegments$CHL_daily_climate <- c(acSegmentsPruned$CHL_daily_climate)

mergedSegments$TKE_surfaceCurrent_5day <- c(acSegmentsPruned$TKE_surfaceCurrent_5day)
mergedSegments$HYCOM_mld_daily <- c(acSegmentsPruned$HYCOM_mld_daily)
mergedSegments$HYCOM_dir_daily <- c(acSegmentsPruned$HYCOM_dir_daily)
mergedSegments$HYCOM_northVel_daily <- c(acSegmentsPruned$HYCOM_northVel_daily)
mergedSegments$HYCOM_eastVel_daily <- c(acSegmentsPruned$HYCOM_eastVel_daily)
mergedSegments$HYCOM_wVel_daily <- c(acSegmentsPruned$HYCOM_wVel_daily)
mergedSegments$Month <- c(acSegmentsPruned$Month)
mergedSegments$SST_8day_climate <- c(acSegmentsPruned$SST_8day_climate)
mergedSegments$SSH_8day_climate <- c(acSegmentsPruned$SSH_8day_climate)
mergedSegments$SST_Monthly_climate <- c(acSegmentsPruned$SST_Monthly_climate)
mergedSegments$SSH_Monthly_climate <- c(acSegmentsPruned$SSH_Monthly_climate)
mergedSegments$EddyDist <- c(acSegmentsPruned$EddyDist)
mergedSegments$Dist_to_Front <- c(acSegmentsPruned$Dist_to_front)
mergedSegments$Weights <- c(rep(1,times = length(acSegmentsPruned$Dist_to_front)))

# mergedSegments <- NULL
# mergedSegments$date <- c(as.Date(visSeg_OnEffort$date_Converted,"%m/%d/%Y"),
#                          as.Date(acSegmentsPruned$date_Converted,"%m/%d/%Y")) #date
# mergedSegments$lat <- c(visSeg_OnEffort$Lat,acSegmentsPruned$Lat)
# mergedSegments$long <- c(visSeg_OnEffort$Long,acSegmentsPruned$Long)
# mergedSegments$ESW <- c(visSeg_OnEffort$ESW,acSegmentsPruned$BW_ESW)
# mergedSegments$SpPresent <- c(visSeg_OnEffort$sp_present,acSegmentsPruned$BW_Present)
# mergedSegments$SpEncounter <- c(visSeg_OnEffort$sp_count,acSegmentsPruned$BW_Encounter)
# mergedSegments$SpEncounter_g0adj <- c(visSeg_OnEffort$sp_count_g0adj,acSegmentsPruned$BW_Encounter)
# mergedSegments$EffectiveArea <- c(visSeg_OnEffort$EffectiveArea,acSegmentsPruned$BW_ESW)
# 
# 
# mergedSegments$Bathymetry <- c(visSeg_OnEffort$Bathymetry,acSegmentsPruned$Bathymetry)
# mergedSegments$SST_daily <- c(visSeg_OnEffort$SST_daily,acSegmentsPruned$SST_daily)
# mergedSegments$SSH_daily <- c(visSeg_OnEffort$SSH_daily,acSegmentsPruned$SSH_daily)
# mergedSegments$CHL_daily_climate <- c(visSeg_OnEffort$CHL_daily_climate,acSegmentsPruned$CHL_daily_climate)
# 
# mergedSegments$TKE_surfaceCurrent_5day <- c(visSeg_OnEffort$TKE_surfaceCurrent_5day,acSegmentsPruned$TKE_surfaceCurrent_5day)
# mergedSegments$HYCOM_mld_daily <- c(visSeg_OnEffort$HYCOM_mld_daily,acSegmentsPruned$HYCOM_mld_daily)
# mergedSegments$HYCOM_dir_daily <- c(visSeg_OnEffort$HYCOM_dir_daily,acSegmentsPruned$HYCOM_dir_daily)
# mergedSegments$HYCOM_northVel_daily <- c(visSeg_OnEffort$HYCOM_northVel_daily,acSegmentsPruned$HYCOM_northVel_daily)
# mergedSegments$HYCOM_eastVel_daily <- c(visSeg_OnEffort$HYCOM_eastVel_daily,acSegmentsPruned$HYCOM_eastVel_daily)
# mergedSegments$HYCOM_wVel_daily <- c(visSeg_OnEffort$HYCOM_wVel_daily,acSegmentsPruned$HYCOM_wVel_daily)
# mergedSegments$Month <- c(visSeg_OnEffort$Month,acSegmentsPruned$Month)
# mergedSegments$SST_8day_climate <- c(visSeg_OnEffort$SST_8day_climate,acSegmentsPruned$SST_8day_climate)
# mergedSegments$SSH_8day_climate <- c(visSeg_OnEffort$SSH_8day_climate,acSegmentsPruned$SSH_8day_climate)
# mergedSegments$SST_Monthly_climate <- c(visSeg_OnEffort$SST_Monthly_climate,acSegmentsPruned$SST_Monthly_climate)
# mergedSegments$SSH_Monthly_climate <- c(visSeg_OnEffort$SSH_Monthly_climate,acSegmentsPruned$SSH_Monthly_climate)
# mergedSegments$EddyDist <- c(visSeg_OnEffort$EddyDist,acSegmentsPruned$EddyDist)
# mergedSegments$Dist_to_Front <- c(visSeg_OnEffort$Dist_to_front,acSegmentsPruned$Dist_to_front)
# 
# mergedSegments$Weights <- c(rep(.1,times = length(visSeg_OnEffort$Dist_to_front)),rep(1,times = length(acSegmentsPruned$Dist_to_front)))

mergedSegments <- as.data.frame(mergedSegments)
covarList<-colnames(mergedSegments[9:length(mergedSegments)])
#covarList<- c("Bathymetry","SST_daily","CHl_8Day","HYCOM_dir_daily","HYCOM_mag_daily","HYCOM_wVel_daily",
#              "SSH_daily","CHL_daily","TKE_surfaceCurrent_5day","HYCOM_mld_daily")
###################################
# Oceanographic variables

# Explore the data, graphical output
cat("Begin exploratory plot generation\n")
# histograms of missing data
percFilled <- plot.missingdata(mergedSegments,covarList) 

# If you decide from the missing data plots that you want to restrict years going forward:
yearListIdx = as.numeric(format(mergedSegments$date,"%Y"))
keepDates.train <- which(yearListIdx != 2012 & yearListIdx >= 2003)
keepDates.test <- which(yearListIdx == 2012)

mergedTrain.set<- mergedSegments[keepDates.train,]

mergedTest.set<- mergedSegments[keepDates.test,]
  

# Cleveland dot plots:
# no transforms
plot.cleveland(mergedTrain.set,covarList,FALSE)


# with transformations
# decided to exclude CHL8day (bad distribution), TKE surface current(outliers, plus redundant),
# SST Monthly climate (Looks the same as 8day climate), SSH Monthly climate (same as 8 day climate),
# bathymetry (not normally distributed for fixed sites)
covarList2 <- c("SST_daily","SSH_daily", "CHL_daily_climate", "HYCOM_mld_daily", "HYCOM_dir_daily",
                "HYCOM_northVel_daily","HYCOM_eastVel_daily", "HYCOM_wVel_daily", "Month",
                "EddyDist","Dist_to_Front")
transformList <- c("none","none","log10","log10","none","none","none","none","none","sqrt","log10")

# restrict covariates again to limited set
mergedTrain.set2<- mergedTrain.set[,covarList2]
mergedTest.set2<- mergedTest.set[,covarList2]

# Identify problematic outliers
outlierList <-which(mergedTrain.set2$Dist_to_Front>300000)
mergedTrain.set2$Dist_to_Front[outlierList] <- NaN 


transformedCovars.train <- transform.covars(mergedTrain.set2,covarList2,transformList)
transformedCovars.test <- transform.covars(mergedTest.set2,covarList2,transformList)

plot.cleveland(transformedCovars.train,colnames(transformedCovars.train),TRUE)

# presence absence histograms
plot.covarDensity(transformedCovars.train,colnames(transformedCovars.train),mergedTrain.set$SpPresent)

# correlation
png(paste(outDir,SP,'_correlations_noTransform.png',sep=''), width = 1000, height = 800)
pairs.panels(transformedCovars.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off()

png(paste(outDir,SP,'_correlations_withTransform.png',sep=''), width = 1000, height = 800)
pairs.panels(transformedCovars.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off() 

cat("Exploratory plots done\n")

###########################
# Run & evaluate models

# # Presence absence
# yBinomial <- mergedTrain.set$SpPresent
# 
# cat("Run full GAM on presence absence data with shrinkage\n")
# 
# presAbsGAMAll <- gam(yBinomial ~ s(SST_daily, bs="ts",k=5) + s(SSH_daily,bs="ts",k=5)
#                     +s(log10_CHL_daily_climate,bs="ts",k=5) + s(log10_HYCOM_mld_daily,bs="ts",k=5) 
#                     +s(HYCOM_northVel_daily,bs="ts",k=5) + s(HYCOM_eastVel_daily,bs="ts",k=5) 
#                     +s(HYCOM_wVel_daily,k=5) + s(EddyDist, bs ="ts",k=5)
#                     +s(log10_Dist_to_Front, bs="ts",k=5)+ s(HYCOM_dir_daily,bs="ts",k=5) ,
#                     method = "GCV.Cp", data = transformedCovars.train, family = binomial(),
#                     offset = log(mergedTrain.set$EffectiveArea),na.action = na.omit)# 
# 
# # Save summary to text file
# sink(paste(outDir,SP,'_GAM_presence_full.txt'))
# summary(presAbsGAMAll)
# gam.check(presAbsGAMAll)
# sink()
# 
# # Calculate and save residuals to text file
# rsd <-residuals.gam(presAbsGAMAll)
# png(paste(outDir,SP,'_residuals_presence_full.png',sep=''), width = 1000, height = 800)
# plot(rsd)
# dev.off() 

# Density
y <- mergedTrain.set$SpEncounter_g0adj
plot(y)
cat("Run full GAM on density data with shrinkage\n")

encounterAll <- gam(y ~ s(SST_daily, bs="ts",k=5) + s(SSH_daily,bs="ts",k=5)
         +s(log10_CHL_daily_climate,bs="ts",k=5) + s(log10_HYCOM_mld_daily,bs="ts",k=5) 
         +s(HYCOM_northVel_daily,bs="ts",k=5)
         +s(EddyDist, bs ="ts",k=5)+ s(Month,bs="ts",k=5)
         +s(log10_Dist_to_Front, bs="ts",k=5)+ s(HYCOM_dir_daily,bs="ts",k=5) ,
         method = "GCV.Cp", data = transformedCovars.train, family = tw(),
         offset = log(mergedTrain.set$EffectiveArea),na.action = na.omit, weights = mergedTrain.set$Weights)# 

# Output summary text to file
sink(paste(outDir,SP,'_GAM_density_full.txt'))
summary(encounterAll)
gam.check(encounterAll)
sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll)
png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 


cat("Run reduced GAM with without unused variables\n")
# look at that output, some covariates have been shrunk down to nothing, so remove them
ctl <- gam.control()
ctl$keepData = TRUE
encounterAll <- gam(y ~ s(SSH_daily, bs="ts",k=5)
                    + s(Month, bs="cc", k=5) 
                    + s(SST_daily, bs="ts", k=5)
                    + s(HYCOM_dir_daily, bs="cc",k=5) ,
                    control = list(keepData=TRUE),
                    method = "GCV.Cp", data = transformedCovars.train, family = tw(),
                    offset = log(mergedTrain.set$EffectiveArea),na.action = na.omit)# 
model <- encounterAll
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, y, NULL, NULL, coordinateSystem, model)
  
plot(encounterAll,pages=1,ylim =c(-2,2))

# Output summary text to file
sink(paste(outDir,SP,'_GAM_density_pruned.txt'))
summary(encounterAll)
gam.check(encounterAll)
sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll)
png(paste(outDir,SP,'_residuals_density_pruned.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 

# save model
save(model, modelMetadata, file = paste(outDir,SP,'_GAM_density_pruned.Rdata',sep=''))

# Predict on test data
yTest <- mergedTest.set$SpEncounter_g0adj
pred <- predict.gam(encounterAll,transformedCovars.test, type = 'response',na.action = na.omit)
plot(yTest,pred)

# Density 
