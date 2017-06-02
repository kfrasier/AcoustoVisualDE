# setup_info_Zc
setwd('E:/NASData/ModelData')

# Species name
SP <- "Ssp" 

#File paths 
acousticSegFile <- "E:/NASData/NASData/ALLSITES_segments_daily.csv" # acoustic input file
acousticDensityFile <-"E:/NASData/MC_GC_DT_binsize000800_Group_density_Cuviers.csv" # acoustic input file"E:/NASData/ALLSITES_binsize000800_Gg_density_jahStart.csv"# 
# The name of the acoustic density file with matched segments
acDensityFile <- paste0('ACDensity_Segments_',SP,'.Rdata')

visDataFile <- "E:/NASData/Sightings_merged.Rdata" # visual sightings
visSegmentsFile <- "E:/NASData/Visual_Segments_v2.csv" # visual segments

# Mapping data
surveyAreaFile <- "E:/NASData/AcoustoVisualDE/surveyAreaOutline.shp"
SPC_vis <- c("Atlantic spotted dolphin", "Striped dolphin","Pantropical spotted dolphin",
             "Spinner dolphin","Stenella sp.","Clymene dolphin")
# Platform Codes (visual only)
PLC <- c("GU","OR")

# Calculate detection functions? This is slow, so if it's already done, you can load trunc dists from file
runDetFuns <- FALSE # can be true or false

# The name of the visual detection function file. 
detFunFile <- paste0("Vis_TruncDist_",SP,"_only.Rdata")# 
# If runDetFuns = TRUE, detFunFile is used to name the R output from the detection function calculation process.
# If runDetFuns = FALSE, detFunFile is used to retrieve the previously caclualated detection functions.

matchACSegs <- FALSE  # set to true if you need to match density estimates with environmental parameters

visG0 <- mean(c(.42,.37)) # for beaked whale from palka 2006 table 5, 2004 survey
              
save(file = "setup_info_Ssp.Rdata", SP,acDensityFile,visDataFile,visSegmentsFile,surveyAreaFile,SPC_vis,PLC,
                                        runDetFuns,detFunFile,matchACSegs,visG0)