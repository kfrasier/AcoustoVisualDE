# SSP models


############################# Acoustic Only Model Fitting #############################
# Run & evaluate models

kVal <- 5
cat("\n Run acoustic only models \n")

yAcOnly_TF <- as.logical(Train_AcOnly.set$Density >0)

cat("Run full binomial GAMs on Acoustic only data with shrinkage\n")#random = list(fac1=~1),
gam_full_AcOnly_TF<-NULL
gam_full_AcOnly_TF$v01 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v02 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1))#
cat("done with model 2: binomial, corAR1, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v03 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v04 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1))#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")
## try with random effects structure
cat("Try adding a random effects structure \n")

gam_full_AcOnly_TF$v05 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 5: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v06 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 6: binomial, corAR1, random effects, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v07 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gam_full_AcOnly_TF$v08 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim'),
                               correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")
# 
# gam_full_AcOnly_TF$v09 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal) 
#                                + s(SSH, bs="re", k=kVal)
#                                + s(log10_CHL, bs="re", k=kVal)
#                                + s(EddyDist, bs="re", k=kVal)
#                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                + s(HYCOM_SALIN_0,bs="re", k=kVal)
#                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                data = transformedCovars_AcOnly.train,
#                                na.action = na.omit,family = nb(),
#                                control = list(opt='optim'),
#                                correlation = corAR1(form=~1|fac1))#
# cat("done with model 9: neg binom, corAR1, EddyDist, 're' spline \n")
# 
# gam_full_AcOnly_TF$v10 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal) 
#                                + s(SSH, bs="re", k=kVal)
#                                + s(log10_CHL, bs="re", k=kVal)
#                                + s(EddyDist, bs="re", k=kVal)
#                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                + s(HYCOM_SALIN_0,bs="re", k=kVal)
#                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                data = transformedCovars_AcOnly.train,
#                                na.action = na.omit,family = binomial(),
#                                control = list(opt='optim'),
#                                correlation = corAR1(form=~1|fac1))#
# cat("done with model 10: binomial, corAR1, EddyDist, 're' spline \n")
# 
# gam_full_AcOnly_TF$v11 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal) 
#                                + s(SSH, bs="re", k=kVal)
#                                + s(log10_CHL, bs="re", k=kVal)
#                                + s(Neg_EddyDist, bs="re", k=kVal)
#                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                + s(HYCOM_SALIN_0,bs="re", k=kVal)
#                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                data = transformedCovars_AcOnly.train,
#                                na.action = na.omit,family = nb(),
#                                control = list(opt='optim'),
#                                correlation = corAR1(form=~1|fac1))#
# cat("done with model 11: neg binom, corAR1, EddyDist, 're' spline \n")
# 
# gam_full_AcOnly_TF$v12 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal) 
#                                + s(SSH, bs="re", k=kVal)
#                                + s(log10_CHL, bs="re", k=kVal)
#                                + s(Neg_EddyDist, bs="re", k=kVal)
#                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                + s(HYCOM_SALIN_0,bs="re", k=kVal)
#                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                data = transformedCovars_AcOnly.train,
#                                na.action = na.omit,family = binomial(),
#                                control = list(opt='optim'),
#                                correlation = corAR1(form=~1|fac1))#

cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")

model_AIC <- NULL
for (iMod in 1:length(gam_full_AcOnly_TF)){
  model_AIC[iMod] = AIC(gam_full_AcOnly_TF[[iMod]]$lme)
  cat('Model ', names(gam_full_AcOnly_TF[iMod]), ' AIC = ', model_AIC[iMod], '\n')
}

best_AcOnly_model <- names(gam_full_AcOnly_TF[which.min(model_AIC)])


cat("\n BEST AC ONLY MODEL IS # ", best_AcOnly_model," \n")
cat("Saving all models to file \n")

save(gam_full_AcOnly_TF, best_AcOnly_model, transformedCovars_AcOnly.train,
     transformedCovars_AcOnly.test,Test_AcOnly.set,Train_AcOnly.set,
     file = paste(outDir,SP,'_AcOnly_binomial_GAMMs_ALL.Rdata',sep=''))

######################### Save best Ac only model ##########################
Ssp_ACOnly_gam_pruned_TF_best <- gamm(yAcOnly_TF ~ s(SST, bs="ts", k=kVal)
                               + s(SSH,bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim', keepData = TRUE),
                               correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
Ssp_ACOnly_gam_pruned_TF_best <- gamm(yAcOnly_TF ~ s(log10_CHL, bs="ts", k=kVal)
                                      + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                      + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                      data = transformedCovars_AcOnly.train,
                                      na.action = na.omit,family = nb(),
                                      control = list(opt='optim', keepData = TRUE),
                                      correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
model <- Ssp_ACOnly_gam_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                                  transformedCovars_AcOnly.train[,c("SST","SSH","log10_CHL",
                                  "HYCOM_SALIN_0","log10_HYCOM_MLD")],
                                  NULL, yAcOnly_TF, NULL, NULL, coordinateSystem, model)
# Then save 'model' to file (.Rdata)
save(model, modelMetadata, file = paste(outDir,SP,'_AcOnly_gamm_pruned_TF_best.Rdata',sep='')) 

# Predict on test data at MC, for comparison...
compAcSet <- which((Test_AcOnly.set$fac1)==21)
pred <- predict.gam(gam_full_AcOnly_TF$v05$gam,transformedCovars_AcOnly.test[compAcSet,],
                    type = 'response',na.action = na.omit)
twoord.plot(Test_AcOnly.set$date[compAcSet],pred,
       Test_AcOnly.set$date[compAcSet],(Test_AcOnly.set$Density[compAcSet]),xlab="Date",type=c("l","b"),
       main="Test Data vs. Model Prediction", ylab = "Model Probability of Occurrence",rylab = "Density From Data")


############################# Visual Only Model Fitting #############################

yVisOnly_TF <- as.logical(Train_VisOnly.set$Density>0)

cat("Run full binomial GAMs on Visual only data with shrinkage\n") # random = list(fac1=~1),
gam_full_VisOnly_TF <- NULL
gam_full_VisOnly_TF$v01 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, Neg_EddyDist, ts spline \n")

gam_full_VisOnly_TF$v02 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 2: binomial, corAR1, Neg_EddyDist, ts spline \n")

gam_full_VisOnly_TF$v03 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline\n")

gam_full_VisOnly_TF$v04 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)                               
                                + s(HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")

## try with random effects structure
cat("Try adding a random effects structure \n")

gam_full_VisOnly_TF$v05 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 5: neg binom, corAR1, random effects, Neg_EddyDist, ts spline \n")

gam_full_VisOnly_TF$v06 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal),
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 6: binomial, corAR1, random effects, Neg_EddyDist, ts spline \n")

gam_full_VisOnly_TF$v07 <- gamm(yVisOnly_TF ~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gam_full_VisOnly_TF$v08 <- gamm(yVisOnly_TF ~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")


## try with random effects smooths
cat("Try random effects smooths instead of thin plate \n")

gam_full_VisOnly_TF$v09 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
                                + s(SSH, bs="re", k=kVal)
                                + s(log10_CHL, bs="re", k=kVal)
                                + s(Neg_EddyDist, bs="re", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                                + s(HYCOM_SALIN_0, bs="re", k=kVal)
                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 9: neg binom, corAR1, Neg_EddyDist, 're' spline \n")

gam_full_VisOnly_TF$v10 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
                                + s(SSH, bs="re", k=kVal)
                                + s(log10_CHL, bs="re", k=kVal)
                                + s(Neg_EddyDist, bs="re", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                                + s(HYCOM_SALIN_0, bs="re", k=kVal)
                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 10: binomial, corAR1, Neg_EddyDist, 're' spline \n")


gam_full_VisOnly_TF$v11 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
                                + s(SSH, bs="re", k=kVal)
                                + s(log10_CHL, bs="re", k=kVal)
                                + s(EddyDist, bs="re", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                                + s(HYCOM_SALIN_0, bs="re", k=kVal)
                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 11: neg binom, corAR1, EddyDist, 're' spline \n")

gam_full_VisOnly_TF$v12 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
                                + s(SSH, bs="re", k=kVal)
                                + s(log10_CHL, bs="re", k=kVal)
                                + s(EddyDist, bs="re", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                                + s(HYCOM_SALIN_0, bs="re", k=kVal)
                                + s(log10_HYCOM_MLD, bs="re", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim'),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")

model_AIC <- NULL
for (iMod in 1:length(gam_full_VisOnly_TF)){
  model_AIC[iMod] = AIC(gam_full_VisOnly_TF[[iMod]]$lme)
  cat('Model ', names(gam_full_VisOnly_TF[[iMod]]), ' AIC = ', model_AIC[iMod], '\n')
}

best_VisOnly_model <- names(gam_full_VisOnly_TF[which.min(model_AIC)])

cat("\n BEST VIS ONLY MODEL IS # ", best_VisOnly_model," \n")
cat("Saving all models to file \n \n")

save(gam_full_VisOnly_TF, best_VisOnly_model, transformedCovars_VisOnly.train,
     transformedCovars_VisOnly.test,Test_VisOnly.set,Train_VisOnly.set,
     file = paste(outDir,SP,'_VisOnly_binomial_GAMMs_ALL.Rdata',sep=''))


#################### Save best vis only model #################################################
Ssp_VisOnly_gamm_pruned_TF <- gamm(yVisOnly_TF~ s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim', keepData = TRUE),
                                correlation = corAR1(form=~1|fac1))

model <- Ssp_VisOnly_gamm_pruned_TF$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                                  transformedCovars_VisOnly.train[,c("HYCOM_SALIN_0","log10_HYCOM_MLD")],
                                  NULL, yVisOnly_TF, NULL, NULL, coordinateSystem, model)


# Then save 'model' to file (.Rdata)
save(model, modelMetadata, file = paste(outDir,SP,'_VisOnly_gamm_pruned_TF_best.Rdata',sep='')) 


########################### Combined Visual and Acoustic Model Fitting ##########################

cat("Run full TF GAMM on combined Visual and Acoustic Data\n")

# Compute weights
myCat <- as.factor(mergedTrain.set$Category)
myWeights <- rep(0,times = length(myCat))
AcTrainSize <- length(which(myCat==2))
VisTrainSize <- length(which(myCat==1))
if (is.null(weight_Ac)|is.null(weight_Vis)){
  # give equal weight to the two datasets
  if (AcTrainSize>=VisTrainSize){
    weight_Vis = round(AcTrainSize/VisTrainSize)
    weight_Ac = 1
  } else {
    weight_Vis = 1
    weight_Ac = round(VisTrainSize/AcTrainSize)
  }
}
myWeights[which(myCat==1)] <- weight_Vis
myWeights[which(myCat==2)] <- weight_Ac

y_TF <- as.logical(mergedTrain.set$Density>0)

cat("Run full binomial GAMs on Visual only data with shrinkage\n")
gam_full_TF <- NULL
gam_full_TF$v01 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train, weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, Neg_EddyDist, ts spline \n")

gam_full_TF$v02 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train, weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100, msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 30)#
cat("done with model 2: binomial, corAR1, Neg_EddyDist, ts spline \n")

gam_full_TF$v03 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline\n")

gam_full_TF$v04 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")

## try with random effects structure
cat("Try adding a random effects structure \n")

gam_full_TF$v05 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 5: neg binom, corAR1, random effects, Neg_EddyDist, ts spline \n")

gam_full_TF$v06 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 50)
cat("done with model 6: binomial, corAR1, random effects, Neg_EddyDist, ts spline \n")

gam_full_TF$v07 <- gamm(y_TF ~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 20)
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gam_full_TF$v08 <- gamm(y_TF ~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 20)
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")


## try with random effects smooths
cat("Try random effects smooths instead of thin plate \n")

gam_full_TF$v09 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
                        + s(SSH, bs="re", k=kVal)
                        + s(log10_CHL, bs="re", k=kVal)
                        + s(Neg_EddyDist, bs="re", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                        + s(HYCOM_SALIN_0, bs="re", k=kVal)
                        + s(log10_HYCOM_MLD, bs="re", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)
cat("done with model 9: neg binom, corAR1, Neg_EddyDist, 're' spline \n")

gam_full_TF$v10 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
                        + s(SSH, bs="re", k=kVal)
                        + s(log10_CHL, bs="re", k=kVal)
                        + s(Neg_EddyDist, bs="re", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                        + s(HYCOM_SALIN_0, bs="re", k=kVal)
                        + s(log10_HYCOM_MLD, bs="re", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100, msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)
cat("done with model 10: binomial, corAR1, Neg_EddyDist, 're' spline \n")

gam_full_TF$v11 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
                        + s(SSH, bs="re", k=kVal)
                        + s(log10_CHL, bs="re", k=kVal)
                        + s(EddyDist, bs="re", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                        + s(HYCOM_SALIN_0, bs="re", k=kVal)
                        + s(log10_HYCOM_MLD, bs="re", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100, msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)
cat("done with model 11: neg binom, corAR1, EddyDist, 're' spline \n")

gam_full_TF$v12 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
                        + s(SSH, bs="re", k=kVal)
                        + s(log10_CHL, bs="re", k=kVal)
                        + s(EddyDist, bs="re", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                        + s(HYCOM_SALIN_0, bs="re", k=kVal)
                        + s(log10_HYCOM_MLD, bs="re", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)
cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")

model_AIC <- NULL
for (iMod in 1:length(gam_full_TF)){
  model_AIC[iMod] = AIC(gam_full_TF[[iMod]]$lme)
  cat('Model ', names(gam_full_TF[iMod]), ' AIC = ', model_AIC[iMod], '\n')
}

best_Combined_model <- names(gam_full_TF[which.min(model_AIC)])

cat("\n BEST COMBINED MODEL IS # ", best_Combined_model," \n")
cat("Saving all models to file \n \n")

save(gam_full_TF, best_Combined_model, transformedCovars.train,
     transformedCovars.test,mergedTest.set,mergedTrain.set,
     file = paste(outDir,SP,'AcousticAndVisual_binomial_GAMMs_ALL.Rdata',sep=''))



#################### Save best full model ################################################i know, but #
Ssp_gamm_pruned_TF_best <- gamm(y_TF~ s(SST, bs="re", k=kVal)
                                + s(log10_CHL, bs="re", k=kVal)
                                + s(EddyDist, bs="re", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
                                + s(HYCOM_SALIN_0, bs="re", k=kVal),
                                data = transformedCovars.train,
                                weights = myWeights,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim',maxIter = 100,msMaxIter=100,keepData = TRUE),
                                correlation = corAR1(form=~1|fac1),niterPQL = 20)

model <- Ssp_gamm_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                                  transformedCovars.train[,c("log10_CHL","HYCOM_SALIN_0",
                                  "log10_HYCOM_MAG_0","EddyDist","SST")],
                                  NULL, y_TF, NULL, NULL, coordinateSystem, model)
# Then save 'model' to file (.Rdata)
save(model, modelMetadata, file = paste(outDir,SP,'_gamm_pruned_TF_best.Rdata',sep='')) 

############# Other Recipes #########

# - To keep data needed for map predictions, add the following cotrol item to gamm call
#       control = list(keepData=TRUE)
# - After running the gam: 
# model <- encounterAll_gam
# coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
# modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, y, NULL, NULL, coordinateSystem, model)
# - Then save 'model' to file (.Rdata)
# save(model, modelMetadata, file = paste(outDir,SP,'_GAM_VisOnly_TF_pruned.Rdata',sep='')) 


# - To output model summary text to file
# sink(paste(outDir,SP,'_GAMM_density_full.txt'))
# summary(encounterAll$gam)
# gam.check(encounterAll$gam)
# sink()


# - To calculate, plot and save residuals to text file
# rsd <-residuals.gam(encounterAll)
# png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
# plot(rsd)
# dev.off() 
# 
# - To plot gam smooths on fewer pages:
# plot(encounterAll$gam,pages=1,ylim =c(-2,2))

