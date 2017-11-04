# spatial prediction: build raster 2009 stacks


#Load in rasters for prediction periods (2009 winter and summer)

predTemplate <-raster('E:/NASData/AcoustoVisualDE/Prediction_template/Predictiontemplate.tif')
log10_CHL_jan_2009 <- log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_jan2009.tif'),
                                     predTemplate))
log10_CHL_july_2009 <- log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_july2009.tif'),
                                      predTemplate))
SST_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_jan_2009.tif'),
                         predTemplate)
SST_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_july_2009.tif'),
                          predTemplate)
SSH_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_jan2009.tif'),
                         predTemplate)
SSH_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_july2009.tif'),
                          predTemplate)
log10_Cayula_jan_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_jan2009.tif'),
           predTemplate)
log10_Cayula_july_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_july2009.tif'),
           predTemplate)
SAL_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_jan2009.tif'),
                         predTemplate)
SAL_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_july2009.tif'),
                          predTemplate)
Eddy_jan_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_Jan_2009_proj.tif'),
                          predTemplate)
Eddy_july_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_July_2009_proj.tif'),
                           predTemplate)
Neg_Eddy_jan_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/proj_negEddyDist_Jan_2009.tif'),
                              predTemplate)
Neg_Eddy_july_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/proj_negEddyDist_Jul_2009.tif'),
                               predTemplate)

Pos_Eddy_jan_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/proj_posEddyDist_Jan_2009.tif'),
                              predTemplate)
Pos_Eddy_july_2009 <- resample(raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/proj_posEddyDist_Jul_2009.tif'),
                               predTemplate)

log10_HYCOM_MAG_0_jan2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_jan2009.tif'),
           predTemplate)
log10_HYCOM_MAG_0_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_july2009.tif'),
           predTemplate)

log10_HYCOM_MLD_jan2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_jan2009.tif'),
           predTemplate)
log10_HYCOM_MLD_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_july2009.tif'),
           predTemplate)

HYCOM_UPVEL_50_jan2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_jan09.tif'),
           predTemplate)
HYCOM_UPVEL_50_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_july09.tif'),
           predTemplate)

jan2009_rasters <- brick(log10_CHL_jan_2009,SST_jan_2009,SSH_jan_2009,log10_Cayula_jan_2009,
                         SAL_jan_2009,Eddy_jan_2009,Neg_Eddy_jan_2009,Pos_Eddy_jan_2009,
                         log10_HYCOM_MAG_0_jan2009,log10_HYCOM_MLD_jan2009,HYCOM_UPVEL_50_jan2009) 

july2009_rasters <- brick(log10_CHL_july_2009,SST_july_2009,SSH_july_2009,log10_Cayula_july_2009,
                          SAL_july_2009,Eddy_july_2009,Neg_Eddy_july_2009,Pos_Eddy_july_2009,
                          log10_HYCOM_MAG_0_july2009,log10_HYCOM_MLD_july2009,HYCOM_UPVEL_50_july2009)      

names(jan2009_rasters)<-c('log10_CHL','SST','SSH','log10_FrontDist_Cayula',
                          'HYCOM_SALIN_0','EddyDist','Neg_EddyDist','Pos_EddyDist',
                          'log10_HYCOM_MAG_0','log10_HYCOM_MLD','HYCOM_UPVEL_50')
names(july2009_rasters)<-c('log10_CHL','SST','SSH','log10_FrontDist_Cayula',
                           'HYCOM_SALIN_0','EddyDist','Neg_EddyDist','Pos_EddyDist',
                           'log10_HYCOM_MAG_0','log10_HYCOM_MLD','HYCOM_UPVEL_50')

save(jan2009_rasters,july2009_rasters,
     file = 'E:/NASData/AcoustoVisualDE/AcoustoVisualDE/2009_prediction_rasters.Rdata')  
