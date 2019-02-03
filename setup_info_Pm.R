# setup_info_Kspp

# Species name
SP <- "Pm" 
SPLong <- 'Sperm whale'

# Set file path with species name
savePath <-(file.path('E:/NASData/ModelData',SP))
# if directory doen't exist, make it.
if (!file.exists(savePath)){
  dir.create(savePath)
}
setwd(savePath)

# Site names included:
siteList= c("MC","GC","DT")#,"DC","MP")

#File paths 
acousticSegFile <- "E:/NASData/ALLSITES_segments_daily_20170918.csv" # acoustic input file
#acousticDensityFile <-"E:/NASData/ModelData/Kspp/MC_GC_DT_DC_MP_binsize011000_Group_density_Kogia.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
acousticDensityFile <-"E:/NASData/ModelData/Pm/MC_GC_DT_binsize_clickDensity_Pm3.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')
pOccurenceFile <- 'E:/NASData/ModelData/Pm/ALLSITES_weeklyPOccurrence_Pm_jahStart.csv' # percent ofdays/week they were present 

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2_20170918.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
SPC_vis <- c("Sperm whale")

# Platform Codes (visual only)
PLC <- c("GU")#"OR"

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- TRUE # can be true or false

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# 
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- TRUE  # set to true if you need to match density estimates with environmental parameters

visG0 <- mean(c(0.57,0.46)) # for Sperm whales from palka 2006 table 5, 2004 survey

# Acoustic truncation distance. Should be the distance within which 95% of detections occur.
AcTruncDist <- 6 # km
Ac_pDet <- 0.50

r_sp <-12 # radius over which acoustic probability of detection applies


save(file = "setup_info_Pm.Rdata", SP,acousticSegFile,acousticDensityFile,visDataFile,SPLong,r_sp,Ac_pDet,
     visSegmentsFile,surveyAreaFile,SPC_vis,PLC,savePath,acDensityFile,
     runDetFuns,detFunFile,matchACSegs,visG0,AcTruncDist,siteList,pOccurenceFile)
