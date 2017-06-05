# best model Zc
Zc_gam_pruned_TF_best <- gamm(y_TF~ s(SST, bs="ts", k=kVal)+s(SSH, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(HYCOM_SALIN_100, bs="ts", k=kVal),
                        data = transformedCovars.train, weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,keepData = TRUE),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20) # Numeric_date-min(Numeric_date)

AIC(Zc_gam_pruned_TF_best$lme)

model <- Zc_gam_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
   SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
   UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                 transformedCovars.train[,c("SST","SSH","EddyDist")],
                 NULL, y_TF, NULL, NULL, coordinateSystem, model)
# Then save 'model' to file (.Rdata)
save(model, modelMetadata, file = paste(outDir,SP,'_gamm_pruned_TF_best.Rdata',sep='')) 


###### Best Acoustic Model ############

Zc_gam_AcOnly_pruned_TF_best <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal) 
                               + s(SSH, bs="re", k=kVal)
                               + s(EddyDist, bs="re", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim',keepData = TRUE),
                               correlation = corAR1(form=~1|fac1))#

AIC(Zc_gam_AcOnly_pruned_TF_best$lme)

Ac_model <- Zc_gam_AcOnly_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
   SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
   UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(Ac_model), "mgcv", 
                                  transformedCovars_AcOnly.train[,c("SST","SSH","EddyDist")],
                                  NULL, yAcOnly_TF, NULL, NULL, coordinateSystem, Ac_model)
# Then save 'model' to file (.Rdata)
save(model, modelMetadata, file = paste(outDir,SP,'_gamm_AcOnly_pruned_TF_best.Rdata',sep='')) 


# ROC curve