# setup_info_Zc
# Species name
SP <- "Zc" 

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
acousticSegFile <- "E:/NASData/ALLSITES_segments_daily.csv" # acoustic input file
acousticDensityFile <-"G:/BW_TPWS/MC_GC_DT_binsize011000_Group_density_Cuviers.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
SPC_vis <- c("Cuvier's beaked whale", "unid. Ziphiid")# "unid. Mesoplodont""Gervais' beaked whale", "Beaked Whale","unid. Mesoplodont"

# Platform Codes (visual only)
PLC <- c("GU","OR")

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- FALSE # can be true or false

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# 
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- TRUE  # set to true if you need to match density estimates with environmental parameters

visG0 <- mean(c(0.27,0.31)) # for beaked whale from palka 2006 table 5, 2004 survey

# Relative datatype weights
weight_Vis <- 1 
weight_Ac <- round(44*(12.6/100)) # in combined model, how should the acoustic data be weighted relative to visual? 

# My logic: If ACOUSTIC estimates are daily bins, 
# and VISUAL estimates are from 10km transect segments at 10 knot survey speed (~18.5km/hr) -> 10/18.5 = 0.54 hr
# then 24/0.54 = 44 (i.e. each ACOUSTIC datapoint is equivalent to 44 VISUAL datapoints)
# BUT, the trucation distance is smaller for acoustic -> surfaceAreaAC ~= 2*2*pi, vs surfaceAreaVis = 10*5*2
# Alternative: If either weight is NULL the two datasets are given equal weights based on dataset size
# ie Ac_weight = number_visSamples/number_acSamples

save(file = "setup_info_Zc.Rdata", SP,acousticSegFile,acousticDensityFile,visDataFile,
     visSegmentsFile,surveyAreaFile,SPC_vis,PLC,savePath,acDensityFile,siteList,
     runDetFuns,detFunFile,matchACSegs,visG0,weight_Ac,weight_Vis)