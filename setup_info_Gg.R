# setup_info_Gg

# Species name
SP <- "Gg" 

# Set file path with species name
savePath <-(file.path('E:/NASData/ModelData',SP))
# if directory doen't exist, make it.
if (!file.exists(savePath)){
  dir.create(savePath)
}
setwd(savePath)

# Site names included:
siteList <-c("MC","GC","DT","DC","MP")

#File paths 
acousticSegFile <- "E:/NASData/ALLSITES_segments_daily_20170918.csv" # acoustic input file
acousticDensityFile <-"E:/NASData/ALLSITES_binsize011000_Gg_density_jahStart.csv"# "E:/NASData/GC_DT_binsize000800_Group_density_Zc.csv" # acoustic input file
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2_20170918.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
SPC_vis <- c("Risso's dolphin")

# Platform Codes (visual only)
PLC <- c("GU")

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- TRUE # can be true or false

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# #Vis_TruncDist_allBW.Rdata" # <- With spares data, you could produce a visual detection function using all beaked whales, 
# then  use that here, to only estimate habitat model for Ziphiid, for instance.
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- TRUE  # set to true if you need to match density estimates with environmental parameters

visG0 <- mean(c(0.77,0.84)) # for Gg  from palka 2006 table 5, 2004 survey
              
# Acoustic truncation distance. Should be the distance within which 95% of detections occur.
AcTruncDist <- 2 # km

save(file = "setup_info_Gg.Rdata", SP,acousticSegFile,acousticDensityFile,visDataFile,
     visSegmentsFile,surveyAreaFile,SPC_vis,PLC,savePath,acDensityFile,siteList,
     runDetFuns,detFunFile,matchACSegs,visG0,AcTruncDist)
