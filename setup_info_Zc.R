# setup_info_Zc

# Species name
SP <- "Zc" 
SPLong <- "Cuvier's beaked whale"
# Set file path with species name
savePath <-(file.path('E:/NASData/ModelData',SP))
# if directory doen't exist, make it.
if (!file.exists(savePath)){
  dir.create(savePath)
}
setwd(savePath)

# Site names included:
siteList= c("MC","GC","DT")

#File paths 
acousticSegFile <- "E:/NASData/ALLSITES_segments_daily_20170918.csv" # acoustic input file
acousticDensityFile <-"E:/NASData/ModelData/Zc/MC_GC_DT_binsize011000_Group_density3_Cuviers.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')
pOccurenceFile <- 'E:/NASData/ModelData/Zc/ALLSITES_weeklyPOccurrence_Zc_jahStart.csv' # percent ofdays/week they were present 

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2_20170918.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
SPC_vis <- c("Cuvier's beaked whale", "unid. Ziphiid")# "unid. Mesoplodont""Gervais' beaked whale", "Beaked Whale","unid. Mesoplodont"

# Platform Codes (visual only)
PLC <- c("GU")#"OR"

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- TRUE # can be true or false

r_sp <-4 # radius over which acoustic probability of detection applies

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# 
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- TRUE  # set to true if you need to match density estimates with environmental parameters

visG0 <- mean(c(0.27,0.31)) # for beaked whale from palka 2006 table 5, 2004 survey

# Acoustic truncation distance. Should be the distance within which 95% of detections occur.
AcTruncDist <- 2.5 # km
Ac_pDet<-0.36
save(file = "setup_info_Zc.Rdata", SP,acousticSegFile,acousticDensityFile,visDataFile,
     visSegmentsFile,surveyAreaFile,SPC_vis,PLC,savePath,acDensityFile,SPLong,r_sp,Ac_pDet,
     runDetFuns,detFunFile,matchACSegs,visG0,AcTruncDist,siteList,pOccurenceFile)
