# make unscaled rasters

predTemplate <-raster('E:/NASData/Eddy/RefRaster.tif')

SST_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_jan_2009.tif'),
                               predTemplate)
SST_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_july_2009.tif'),
                                predTemplate)

SSH_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_jan2009.tif'),
                               predTemplate)
SSH_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_july2009.tif'),
                                predTemplate)

log10_CHL_jan_2009 <- log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_jan2009.tif'),
                                           predTemplate))
log10_CHL_july_2009 <- log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_july2009.tif'),
                                            predTemplate))

log10_HYCOM_MLD_jan2009 <- resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_jan2009.tif'),
           predTemplate)
log10_HYCOM_MLD_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_july2009.tif'),
           predTemplate)

SAL0_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_jan2009.tif'),
                                predTemplate)
                       
SAL0_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_july2009.tif'),
                                 predTemplate)
SAL100_jan_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal100_jan2009_albers.tif'),
                                  predTemplate)
SAL100_july_2009 <- resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal100_july2009_albers.tif'),
                                   predTemplate)
log10_HYCOM_MAG_0_jan2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_jan2009.tif'),
           predTemplate)

log10_HYCOM_MAG_0_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_july2009.tif'),
           predTemplate)

HYCOM_UPVEL_jan2009 <-  
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_jan09.tif'),
           predTemplate)

HYCOM_UPVEL_july2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_july09.tif'),
           predTemplate)

log10_Cayula_jan_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_jan2009.tif'),
                 predTemplate)

log10_Cayula_july_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_july2009.tif'),
                 predTemplate)

Eddy_jan_2009 <- resample(
  raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_Jan_2009_proj.tif'),
  predTemplate)

Eddy_july_2009 <- resample(raster(
  'E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_July_2009_proj.tif'),
  predTemplate)

Neg_Eddy_jan_2009 <- resample(raster(
  'E:/NASData/AcoustoVisualDE/EddyDist/DistNegMSLA_eddy_polarities_2009015.img')/1000,
  predTemplate)

Neg_Eddy_july_2009 <- resample(raster(
  'E:/NASData/AcoustoVisualDE/EddyDist/Neg_EddyDist_july2009.tif')/1000,
  predTemplate)

Pos_Eddy_jan_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/EddyDist/PosEddyDist_Jan_2009.tif'),
           predTemplate)

Pos_Eddy_july_2009 <- 
  resample(raster('E:/NASData/AcoustoVisualDE/EddyDist/PosEddyDist_Jul_2009.tif'),
           predTemplate)


jan2009_rasters <- brick(log10_CHL_jan_2009,SST_jan_2009,SSH_jan_2009,log10_Cayula_jan_2009,
                         SAL0_jan_2009,SAL100_jan_2009, 
                         Eddy_jan_2009,Neg_Eddy_jan_2009,Pos_Eddy_jan_2009,
                         log10_HYCOM_MAG_0_jan2009,log10_HYCOM_MLD_jan2009,HYCOM_UPVEL_jan2009) 

july2009_rasters <- brick(log10_CHL_july_2009,SST_july_2009,SSH_july_2009,log10_Cayula_july_2009,
                          SAL0_july_2009,SAL100_july_2009,
                          Eddy_july_2009,Neg_Eddy_july_2009,Pos_Eddy_july_2009,
                          log10_HYCOM_MAG_0_july2009,log10_HYCOM_MLD_july2009,HYCOM_UPVEL_july2009)      

names(jan2009_rasters)<-c('log10_CHL','SST','SSH','log10_FrontDist_Cayula',
                          'HYCOM_SALIN_0','HYCOM_SALIN_100',
                          'EddyDist','Neg_EddyDist',"Pos_EddyDist",
                          'log10_HYCOM_MAG_0','log10_HYCOM_MLD','HYCOM_UPVEL_50')
names(july2009_rasters)<-c('log10_CHL','SST','SSH','log10_FrontDist_Cayula',
                           'HYCOM_SALIN_0','HYCOM_SALIN_100',
                           'EddyDist','Neg_EddyDist',"Pos_EddyDist",
                           'log10_HYCOM_MAG_0','log10_HYCOM_MLD','HYCOM_UPVEL_50')

save(jan2009_rasters,july2009_rasters,
     file = 'E:/NASData/AcoustoVisualDE/AcoustoVisualDE/2009_prediction_rasters_UNscaled.Rdata')  