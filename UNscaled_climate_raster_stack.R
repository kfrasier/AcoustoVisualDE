# UN_scaled_climate_raster_stack 
  
predTemplate <-raster('E:/NASData/Eddy/RefRaster.tif')
#map_proj <- crs(predTemplate)
  
# folders
eddy_ClimateDir <- 'E:/NASData/Eddy/JPL_ManualEddyDist/Climatologies/Projected'
SST_ClimateDir <- 'E:/NASData/SST_monthly_climatology_20030101to20151231/GLOB/JPL_OUROCEAN/G1SST/analysed_sst/Monthly_Climatology/Projected'
SSH_ClimateDir <- 'E:/NASData/SSH_monthly_climatology_20030101to20151231'
CHL_ClimateDir <- 'E:/NASData/CHL_monthly_climatology_20030101to20151231/aqua/monthly/4km/CHL_chlor_a/Monthly_Climatology/Projected'
HYCOM_MAG_ClimateDir <- 'E:/NASData/HYCOM_monthly_climatology_20030101to20151231/mag/Monthly_Climatology/Depth_0000m/Projected'
HYCOM_SALIN_ClimateDir <- 'E:/NASData/HYCOM_monthly_climatology_20030101to20151231/salinity/Monthly_Climatology/Depth_0000m/Projected'
HYCOM_MLD_ClimateDir <- 'E:/NASData/HYCOM_monthly_climatology_20030101to20151231/mld/Monthly_Climatology/Projected'
HYCOM_UPVEL_ClimateDir <-'E:/NASData/HYCOM_monthly_climatology_20030101to20151231/w_velocity/Monthly_Climatology/Depth_0050m/Projected'

monthNum <-c('01','02','03','04','05','06','07','08','09','10','11','12')
monthStr<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
raster_set <- vector('list',length = 12)

for (iM in 1:length(monthStr)){
  SST_climate <- resample(raster(paste0(SST_ClimateDir,'/',
                       'proj_JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST-analysed_sst-month',
                       monthNum[iM],'-mean.img')),predTemplate)
  
  SSH_climate <- resample(raster(paste0(SSH_ClimateDir,'/','proj_SSH_',
                       monthStr[iM],'_climatology.tif')),predTemplate)
  
  log10_CHL_climate <- log10(resample(raster(paste0(CHL_ClimateDir,'/',
                       'proj_aqua_CHL_chlor_a_month',
                       monthNum[iM],'_mean.img')),predTemplate))
  
  log10_HYCOM_MLD_climate <- log10(resample(raster(paste0(HYCOM_MLD_ClimateDir,'/',
                             'proj_mld_month',monthNum[iM],'_mean.img')),predTemplate))
  
  HYCOM_SALIN_climate <- resample(raster(paste0(HYCOM_SALIN_ClimateDir,'/',
                         'proj_salinity_0000m_month',monthNum[iM],'_mean.img')),predTemplate)
  
  log10_HYCOM_MAG_climate <- log10(resample(raster(paste0(HYCOM_MAG_ClimateDir,'/',
                              'proj_mag_0000m_month', monthNum[iM],'_mean.img')),predTemplate))
  
  HYCOM_UPVEL_climate <- resample(raster(paste0(HYCOM_UPVEL_ClimateDir,'/',
                         'proj_w_velocity_0050m_month', monthNum[iM],'_mean.img')),predTemplate)
  
  eddy_climate <- resample(raster(paste0(eddy_ClimateDir,'/','proj_eddyDist_',
                  monthStr[iM],'_climatology.tif')),predTemplate)
  
  pos_eddy_climate <- resample(raster(paste0(eddy_ClimateDir,'/','proj_posEddyDist_',
                      monthStr[iM],'_climatology.tif')),predTemplate)
  
  neg_eddy_climate <- resample(raster(paste0(eddy_ClimateDir,'/','proj_negEddyDist_',
                      monthStr[iM],'_climatology.tif')),predTemplate)
  
  raster_set[[iM]] <- brick(SST_climate,SSH_climate,log10_CHL_climate,eddy_climate,
                            pos_eddy_climate,neg_eddy_climate,log10_HYCOM_MAG_climate,
                            HYCOM_UPVEL_climate,HYCOM_SALIN_climate,log10_HYCOM_MLD_climate) 
  
  names(raster_set[[iM]])<-c('SST','SSH','log10_CHL','EddyDist',
                             'Pos_EddyDist','Neg_EddyDist','log10_HYCOM_MAG_0',
                             'HYCOM_UPVEL_50','HYCOM_SALIN_0','log10_HYCOM_MLD')  
  cat(paste0("done with month ", iM, "\n", collapse = ""))
  
}

save(raster_set,monthStr, monthNum,
     file = 'E:/NASData/AcoustoVisualDE/AcoustoVisualDE/climatology_rasters_UNscaled.Rdata')  

