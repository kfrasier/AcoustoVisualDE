# make SCALED rasters

### this code only needs to be run to regenerate the raster stack.


predTemplate <-raster('E:/NASData/Eddy/RefRaster.tif')

SST_jan_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_jan_2009.tif'),
                               predTemplate),center = covars_Joint_min.train["SST"],
                      scale = covars_Joint_max.train["SST"]-covars_Joint_min.train["SST"])
SST_july_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/SST/SST_july_2009.tif'),
                                predTemplate),center = covars_Joint_min.train["SST"],
                       scale = covars_Joint_max.train["SST"]-covars_Joint_min.train["SST"])

SSH_jan_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_jan2009.tif'),
                               predTemplate),center = covars_Joint_min.train["SSH"],
                      scale = covars_Joint_max.train["SSH"]-covars_Joint_min.train["SSH"])
SSH_july_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/SSH/SSH_july2009.tif'),
                                predTemplate),center = covars_Joint_min.train["SSH"],
                       scale = covars_Joint_max.train["SSH"]-covars_Joint_min.train["SSH"])

log10_CHL_jan_2009 <- scale(log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_jan2009.tif'),
                                           predTemplate)),center = covars_Joint_min.train["log10_CHL"],scale = 
                              covars_Joint_max.train["log10_CHL"]-covars_Joint_min.train["log10_CHL"])
log10_CHL_july_2009 <- scale(log10(resample(raster('E:/NASData/AcoustoVisualDE/CHL/CHL_july2009.tif'),
                                            predTemplate)),center = covars_Joint_min.train["log10_CHL"],scale = 
                               covars_Joint_max.train["log10_CHL"]-covars_Joint_min.train["log10_CHL"])

log10_HYCOM_MLD_jan2009 <- scale( 
  resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_jan2009.tif'),
           predTemplate),center = covars_Joint_min.train["log10_HYCOM_MLD"],
  scale = covars_Joint_max.train["log10_HYCOM_MLD"]-covars_Joint_min.train["log10_HYCOM_MLD"])
log10_HYCOM_MLD_july2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/MLD/log10_mld_july2009.tif'),
           predTemplate),center = covars_Joint_min.train["log10_HYCOM_MLD"],
  scale = covars_Joint_max.train["log10_HYCOM_MLD"]-covars_Joint_min.train["log10_HYCOM_MLD"])

SAL0_jan_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_jan2009.tif'),
                                predTemplate),center = covars_Joint_min.train["HYCOM_SALIN_0"],
                       scale = covars_Joint_max.train["HYCOM_SALIN_0"]-covars_Joint_min.train["HYCOM_SALIN_0"])
SAL0_july_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal0_july2009.tif'),
                                 predTemplate),center = covars_Joint_min.train["HYCOM_SALIN_0"],
                        scale = covars_Joint_max.train["HYCOM_SALIN_0"]-covars_Joint_min.train["HYCOM_SALIN_0"])

SAL100_jan_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal100_jan2009_albers.tif'),
                                  predTemplate),center = covars_Joint_min.train["HYCOM_SALIN_100"],
                         scale = covars_Joint_max.train["HYCOM_SALIN_100"]-covars_Joint_min.train["HYCOM_SALIN_100"])
SAL100_july_2009 <- scale(resample(raster('E:/NASData/AcoustoVisualDE/Salinity/sal100_july2009_albers.tif'),
                                   predTemplate),center = covars_Joint_min.train["HYCOM_SALIN_100"],
                          scale = covars_Joint_max.train["HYCOM_SALIN_100"]-covars_Joint_min.train["HYCOM_SALIN_100"])

log10_HYCOM_MAG_0_jan2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_jan2009.tif'),
           predTemplate),center = covars_Joint_min.train['log10_HYCOM_MAG_0'],
  scale = covars_Joint_max.train['log10_HYCOM_MAG_0']-covars_Joint_min.train['log10_HYCOM_MAG_0'])
log10_HYCOM_MAG_0_july2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_Mag/log10_mag_july2009.tif'),
           predTemplate),center = covars_Joint_min.train['log10_HYCOM_MAG_0'],
  scale = covars_Joint_max.train['log10_HYCOM_MAG_0']-covars_Joint_min.train['log10_HYCOM_MAG_0'])

HYCOM_UPVEL_jan2009 <- scale( 
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_jan09.tif'),
           predTemplate),center = covars_Joint_min.train['HYCOM_UPVEL_50'],
  scale = covars_Joint_max.train['HYCOM_UPVEL_50']-covars_Joint_min.train['HYCOM_UPVEL_50'])
HYCOM_UPVEL_july2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/HYCOM_UpVel/proj_HYCOM_upvel_50_july09.tif'),
           predTemplate),center = covars_Joint_min.train['HYCOM_UPVEL_50'],
  scale = covars_Joint_max.train['HYCOM_UPVEL_50']-covars_Joint_min.train['HYCOM_UPVEL_50'])

log10_Cayula_jan_2009 <- 
  scale(resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_jan2009.tif'),
                 predTemplate),center = covars_Joint_min.train["log10_FrontDist_Cayula"],
        scale = covars_Joint_max.train["log10_FrontDist_Cayula"]-covars_Joint_min.train["log10_FrontDist_Cayula"])
log10_Cayula_july_2009 <- 
  scale(resample(raster('E:/NASData/AcoustoVisualDE/Cayula/log10_CayulaFront_july2009.tif'),
                 predTemplate),center = covars_Joint_min.train["log10_FrontDist_Cayula"],
        scale = covars_Joint_max.train["log10_FrontDist_Cayula"]-covars_Joint_min.train["log10_FrontDist_Cayula"])

Eddy_jan_2009 <- scale(resample(
  raster('E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_Jan_2009_proj.tif'),
  predTemplate),center = covars_Joint_min.train["EddyDist"],
  scale = covars_Joint_max.train["EddyDist"]-covars_Joint_min.train["EddyDist"])
Eddy_july_2009 <- scale(resample(raster(
  'E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected/eddyDist_July_2009_proj.tif'),
  predTemplate),center = covars_Joint_min.train["EddyDist"],
  scale = covars_Joint_max.train["EddyDist"]-covars_Joint_min.train["EddyDist"])

Neg_Eddy_jan_2009 <- scale(resample(raster(
  'E:/NASData/AcoustoVisualDE/EddyDist/DistNegMSLA_eddy_polarities_2009015.img')/1000,
  predTemplate),center = covars_Joint_min.train["Neg_EddyDist"],
  scale = covars_Joint_max.train["Neg_EddyDist"]-covars_Joint_min.train["Neg_EddyDist"])
Neg_Eddy_july_2009 <- scale(resample(raster(
  'E:/NASData/AcoustoVisualDE/EddyDist/Neg_EddyDist_july2009.tif')/1000,
  predTemplate),center = covars_Joint_min.train["Neg_EddyDist"],
  scale = covars_Joint_max.train["Neg_EddyDist"]-covars_Joint_min.train["Neg_EddyDist"])

Pos_Eddy_jan_2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/EddyDist/PosEddyDist_Jan_2009.tif'),
           predTemplate),center = covars_Joint_min.train["Pos_EddyDist"],
  scale = covars_Joint_max.train["Pos_EddyDist"]-covars_Joint_min.train["Pos_EddyDist"])
Pos_Eddy_july_2009 <- scale(
  resample(raster('E:/NASData/AcoustoVisualDE/EddyDist/PosEddyDist_Jul_2009.tif'),
           predTemplate),center = covars_Joint_min.train["Pos_EddyDist"],
  scale = covars_Joint_max.train["Pos_EddyDist"]-covars_Joint_min.train["Pos_EddyDist"])


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
     file = 'E:/NASData/AcoustoVisualDE/AcoustoVisualDE/2009_prediction_rasters_scaled.Rdata')  
