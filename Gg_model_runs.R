# Gg models
library(mgcv)
library(MASS)
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/GetModelMetadata.R')

outDir <- file.path("E:/NASData/ModelData",SP,"/")

load('E:/NASData/ModelData/Gg/setup_info_Gg.Rdata')
weeklyOccurrenceFile = paste0(outDir,'ALLSITES_weeklyPOccurrence_Gg_jahStart.csv')
pOccur <- read.csv(weeklyOccurrenceFile, header = TRUE,na.strings=c(""," ","NA","NaN"))

############################# Acoustic Only Model Fitting #############################
# Run & evaluate models

kVal <- 5
cat("\n Run acoustic only models \n")

yAcOnly_TF <- as.logical(Train_AcOnly.set$Density >0)

cat("Run full binomial GAMs on Acoustic only data with shrinkage\n")#random = list(fac1=~1),
gamm_full_AcOnly_TF<-NULL
gamm_full_AcOnly_TF$v01 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v02 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1))#
cat("done with model 2: binomial, corAR1, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v03 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v04 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1))#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")
## try with random effects structure
cat("Try adding a random effects structure \n")

gamm_full_AcOnly_TF$v05 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 5: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v06 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 6: binomial, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v07 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = nb(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v08 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                               + s(SSH, bs="ts", k=kVal)
                               + s(log10_CHL, bs="ts", k=kVal)
                               + s(Neg_EddyDist, bs="ts", k=kVal)
                               + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                               + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                               + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                               + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                               data = transformedCovars_AcOnly.train,
                               na.action = na.omit,family = binomial(),
                               control = list(opt='optim', niterEM=0),
                               correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")

# gamm_full_AcOnly_TF$v09 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal)
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
# gamm_full_AcOnly_TF$v10 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal)
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
# gamm_full_AcOnly_TF$v11 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal)
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
# gamm_full_AcOnly_TF$v12 <- gamm(yAcOnly_TF~ s(SST, bs="re", k=kVal)
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

# cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")

gamm_full_AcOnly_TF$v13 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_AcOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 13: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v14 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_AcOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 14: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v15 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_AcOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 15: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_AcOnly_TF$v16 <- gamm(yAcOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_AcOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 16: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")

model_AIC <- NULL
for (iMod in 1:length(gamm_full_AcOnly_TF)){
  model_AIC[iMod] = AIC(gamm_full_AcOnly_TF[[iMod]]$lme)
  cat('Model ', names(gamm_full_AcOnly_TF[iMod]), ' AIC = ', model_AIC[iMod], '\n')
}

best_AcOnly_model <- names(gamm_full_AcOnly_TF[which.min(model_AIC)])


cat("\n BEST AC ONLY MODEL IS # ", best_AcOnly_model," \n")
cat("Saving all models to file \n")

save(gamm_full_AcOnly_TF, best_AcOnly_model, transformedCovars_AcOnly.train,
     transformedCovars_AcOnly.test,Test_AcOnly.set,Train_AcOnly.set,
     file = paste(outDir,SP,'_AcOnly_binomial_GAMMs_ALL.Rdata',sep=''))

######################### Save best Ac only model ##########################
ACOnly_gamm_pruned_TF_best <- gamm(yAcOnly_TF ~ s(SST, bs="ts", k=kVal)
                                      + s(SSH, bs="ts", k=kVal)
                                      + s(EddyDist, bs="ts", k=kVal)
                                      + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                      + s(HYCOM_SALIN_0,bs="ts", k=kVal),
                                      data = transformedCovars_AcOnly.train,
                                      na.action = na.omit,family = nb(),
                                      control = list(opt='optim',niterEM=0,keepData = TRUE),
                                      correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)

model <- ACOnly_gamm_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                                  transformedCovars_AcOnly.train[,c("SST","SSH","log10_CHL",
                                  "HYCOM_SALIN_0","EddyDist","log10_FrontDist_Cayula")],
                                  NULL, yAcOnly_TF, NULL, NULL, coordinateSystem, model)
# Then save 'model' to file (.Rdata)
save(model, ACOnly_gamm_pruned_TF_best,
     modelMetadata, file = paste(outDir,SP,'_AcOnly_gamm_pruned_TF_best.Rdata',sep='')) 


# Output best acoustic only model summary text to file
sink(paste0(outDir,SP,'_AcOnly_GAMM_TF_pruned_best_summary.txt'))
summary(model)
sink()

# Plot model smooths
png(paste0(outDir,SP,'_AcOnly_BestModel_smooths.png',sep=''),width = 2000,height = 1600,res=300)
plot(model,pages=1,cex.lab = 1.3,cex.axis = 1.1,scale=0)
dev.off()

# Predict on acoustic test data, using acoustic only model for comparison...
compAcSet_MC <- which((Test_AcOnly.set$fac1)==5)
compAcSet_GC <- which((Test_AcOnly.set$fac1)==10)
compAcSet_DT <- which((Test_AcOnly.set$fac1)==15 |(Test_AcOnly.set$fac1)==16)
compAcSet_DC <- which((Test_AcOnly.set$fac1)==21 |(Test_AcOnly.set$fac1)==22)
compAcSet_MP <- which((Test_AcOnly.set$fac1)==26)

dateTicks = as.POSIXct(c("2013-01-01 GMT","2013-04-01 GMT","2013-07-01 GMT","2013-10-01 GMT","2014-01-01 GMT"))
dateLabels = c("Jan. 2013","Apr. 2013","Jul. 2013","Oct 2013","Jan. 2014")
predAcOnly_MC <- predict.gam(ACOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MC,],
                       type = 'response',na.action = na.pass)
predAcOnly_GC <- predict.gam(ACOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_GC,],
                       type = 'response',na.action = na.pass)
predAcOnly_DT <- predict.gam(ACOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DT,],
                       type = 'response',na.action = na.pass)
predAcOnly_DC <- predict.gam(ACOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DC,],
                       type = 'response',na.action = na.pass)
predAcOnly_MP <- predict.gam(ACOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MP,],
                       type = 'response',na.action = na.pass)
occurIdx = which(as.POSIXct(pOccur[,1])>='2013-01-01' &as.POSIXct(pOccur[,1])<'2014-01-01')

png(paste0(outDir,SP,'_AcOnly_Prediction_timeseries.png',sep=''),width = 2000,height = 3000,res=300)
par(mfrow = c(5,1))
twoord.plot(Test_AcOnly.set$date[compAcSet_MC],predAcOnly_MC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,2],type=c("l","l"),mar = c(1, 5, 1, 5),
            main="Test Data vs. Model Prediction",xtickpos = dateTicks,xticklab = rep("",5),axislab.cex=.8)
twoord.plot(Test_AcOnly.set$date[compAcSet_GC],predAcOnly_GC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,3],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DT],predAcOnly_DT,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,4],type=c("l","l"),
            ylab = "Model Probability of Occurrence",mar = c(1, 5, 0, 5),
            rylab = "Weekly Probability of Occurrence from Data",xtickpos = dateTicks,axislab.cex=.8,
            xticklab = rep("",5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DC],predAcOnly_DC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,5],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_MP],predAcOnly_MP,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,6],type=c("l","l"),mar = c(3, 5, 0, 5),
            xlab = "Date",xtickpos = dateTicks, xticklab = dateLabels,axislab.cex=.8)

#axis(1,at=dateTicks, labels=dateLabels)
dev.off()

############################# Visual Only Model Fitting #############################

yVisOnly_TF <- as.logical(Train_VisOnly.set$Density>0)

cat("Run full binomial GAMMs on Visual only data with shrinkage\n") # random = list(fac1=~1),
gamm_full_VisOnly_TF <- NULL
gamm_full_VisOnly_TF$v01 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1)) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, Neg_EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v02 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 2: binomial, corAR1, Neg_EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v03 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline\n")

gamm_full_VisOnly_TF$v04 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)    
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)                           
                                + s(HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")

## try with random effects structure
cat("Try adding a random effects structure \n")

gamm_full_VisOnly_TF$v05 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 5: neg binom, corAR1, random effects, Neg_EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v06 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 6: binomial, corAR1, random effects, Neg_EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v07 <- gamm(yVisOnly_TF ~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = nb(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v08 <- gamm(yVisOnly_TF ~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = binomial(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1),random=list(fac1=~1))#
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")


# ## try with random effects smooths
# cat("Try random effects smooths instead of thin plate \n")
# 
# gamm_full_VisOnly_TF$v09 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
#                                 + s(SSH, bs="re", k=kVal)
#                                 + s(log10_CHL, bs="re", k=kVal)
#                                 + s(Neg_EddyDist, bs="re", k=kVal)
#                                 + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                 + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                                 + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                 + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                 data = transformedCovars_VisOnly.train,
#                                 na.action = na.omit,family = nb(),
#                                 control = list(opt='optim',niterEM=0),
#                                 correlation = corAR1(form=~1|fac1))#
# cat("done with model 9: neg binom, corAR1, Neg_EddyDist, 're' spline \n")
# 
# gamm_full_VisOnly_TF$v10 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
#                                 + s(SSH, bs="re", k=kVal)
#                                 + s(log10_CHL, bs="re", k=kVal)
#                                 + s(Neg_EddyDist, bs="re", k=kVal)
#                                 + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                 + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                                 + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                 + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                 data = transformedCovars_VisOnly.train,
#                                 na.action = na.omit,family = binomial(),
#                                 control = list(opt='optim',niterEM=0),
#                                 correlation = corAR1(form=~1|fac1))#
# cat("done with model 10: binomial, corAR1, Neg_EddyDist, 're' spline \n")
# 
# 
# gamm_full_VisOnly_TF$v11 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
#                                 + s(SSH, bs="re", k=kVal)
#                                 + s(log10_CHL, bs="re", k=kVal)
#                                 + s(EddyDist, bs="re", k=kVal)
#                                 + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                 + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                                 + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                 + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                 data = transformedCovars_VisOnly.train,
#                                 na.action = na.omit,family = nb(),
#                                 control = list(opt='optim',niterEM=0),
#                                 correlation = corAR1(form=~1|fac1))#
# cat("done with model 11: neg binom, corAR1, EddyDist, 're' spline \n")
# 
# gamm_full_VisOnly_TF$v12 <- gamm(yVisOnly_TF~ s(SST, bs="re", k=kVal) 
#                                 + s(SSH, bs="re", k=kVal)
#                                 + s(log10_CHL, bs="re", k=kVal)
#                                 + s(EddyDist, bs="re", k=kVal)
#                                 + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                                 + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                                 + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                                 + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                                 data = transformedCovars_VisOnly.train,
#                                 na.action = na.omit,family = binomial(),
#                                 control = list(opt='optim',niterEM=0),
#                                 correlation = corAR1(form=~1|fac1))#
# cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")

gamm_full_VisOnly_TF$v13 <- gamm(yVisOnly_TF ~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 13: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v14 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 14: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v15 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim', niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 15: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_VisOnly_TF$v16 <- gamm(yVisOnly_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                data = transformedCovars_VisOnly.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim',niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 16: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")
model_AIC <- NULL
for (iMod in 1:length(gamm_full_VisOnly_TF)){
  model_AIC[iMod] = AIC(gamm_full_VisOnly_TF[[iMod]]$lme)
  cat('Model ', names(gamm_full_VisOnly_TF[iMod]), ' AIC = ', model_AIC[iMod], '\n')
}

best_VisOnly_model <- names(gamm_full_VisOnly_TF[which.min(model_AIC)])

cat("\n BEST VIS ONLY MODEL IS # ", best_VisOnly_model," \n")
cat("Saving all models to file \n \n")

save(gamm_full_VisOnly_TF, best_VisOnly_model, transformedCovars_VisOnly.train,
     transformedCovars_VisOnly.test,Test_VisOnly.set,Train_VisOnly.set,
     file = paste(outDir,SP,'_VisOnly_binomial_GAMMs_ALL.Rdata',sep=''))


#################### Save best vis only model #################################################
VisOnly_gamm_pruned_TF_best <- gamm(yVisOnly_TF~ s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                   data = transformedCovars_VisOnly.train,
                                   na.action = na.omit,family = nb(),
                                   control = list(opt='optim', niterEM=0,keepData = TRUE),
                                   correlation = corAR1(form=~1|fac1),random=list(fac1=~1))

model <- VisOnly_gamm_pruned_TF$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                 transformedCovars_VisOnly.train[,c("SST","log10_HYCOM_MLD")],
                 NULL, yVisOnly_TF, NULL, NULL, coordinateSystem, model)

# Then save 'model' to file (.Rdata)
save(model, modelMetadata, VisOnly_gamm_pruned_TF,
     file = paste(outDir,SP,'_VisOnly_gamm_pruned_TF_best.Rdata',sep='')) 

# Plot model smooths
png(paste0(outDir,SP,'_VisOnly_BestModel_smooths.png',sep=''),width = 1000,height = 1000,res=300)
plot(model,pages=1)
dev.off()

# Output best acoustic only model summary text to file
sink(paste0(outDir,SP,'_VisOnly_GAMM_TF_pruned_best_summary.txt'))
summary(model)
sink()

# Temporal prediction
# Predict on acoustic test data using visual only model, for comparison...
predVisOnly_MC <- predict.gam(VisOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MC,],
                       type = 'response',na.action = na.pass)
predVisOnly_GC <- predict.gam(VisOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_GC,],
                       type = 'response',na.action = na.pass)
predVisOnly_DT <- predict.gam(VisOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DT,],
                       type = 'response',na.action = na.pass)
predVisOnly_DC <- predict.gam(VisOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DC,],
                       type = 'response',na.action = na.pass)
predVisOnly_MP <- predict.gam(VisOnly_gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MP,],
                       type = 'response',na.action = na.pass)
occurIdx = which(as.POSIXct(pOccur[,1])>='2013-01-01' &as.POSIXct(pOccur[,1])<'2014-01-01')

png(paste0(outDir,SP,'_VisOnly_Prediction_timeseries.png',sep=''),width = 2000,height = 3000,res=300)
par(mfrow = c(5,1))
twoord.plot(Test_AcOnly.set$date[compAcSet_MC],predVisOnly_MC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,2],type=c("l","l"),mar = c(1, 5, 1, 5),
            main="Test Data vs. Model Prediction",xtickpos = dateTicks,xticklab = rep("",5),axislab.cex=.8)
twoord.plot(Test_AcOnly.set$date[compAcSet_GC],predVisOnly_GC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,3],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DT],predVisOnly_DT,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,4],type=c("l","l"),
            ylab = "Model Probability of Occurrence",mar = c(1, 5, 0, 5),
            rylab = "Weekly Probability of Occurrence from Data",xtickpos = dateTicks,axislab.cex=.8,
            xticklab = rep("",5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DC],predVisOnly_DC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,5],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_MP],predVisOnly_MP,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,6],type=c("l","l"),mar = c(3, 5, 0, 5),
            xlab = "Date",xtickpos = dateTicks, xticklab = dateLabels,axislab.cex=.8)

#axis(1,at=dateTicks, labels=dateLabels)
dev.off()
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

cat("Run full binomial GAMMs on Visual only data with shrinkage\n")
gamm_full_TF <- NULL
gamm_full_TF$v01 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train, weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter = 100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20) # Numeric_date-min(Numeric_date)
cat("done with model 1: neg binom, corAR1, Neg_EddyDist, ts spline \n")

gamm_full_TF$v02 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train, weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100, msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 2: binomial, corAR1, Neg_EddyDist, ts spline \n")

gamm_full_TF$v03 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 3: neg binom, corAR1, EddyDist, ts spline\n")

gamm_full_TF$v04 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 4: binomial, corAR1, EddyDist, ts spline \n")

## try with random effects structure
cat("Try adding a random effects structure \n")

gamm_full_TF$v05 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)#
cat("done with model 5: neg binom, corAR1, random effects, Neg_EddyDist, ts spline \n")

gamm_full_TF$v06 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts",k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 50)
cat("done with model 6: binomial, corAR1, random effects, Neg_EddyDist, ts spline \n")

gamm_full_TF$v07 <- gamm(y_TF ~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 20)
cat("done with model 7: neg binom, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_TF$v08 <- gamm(y_TF ~ s(SST, bs="ts", k=kVal) 
                        + s(SSH, bs="ts", k=kVal)
                        + s(log10_CHL, bs="ts", k=kVal)
                        + s(EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                        + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = binomial(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                        correlation = corAR1(form=~1|fac1),
                        random=list(fac1=~1),niterPQL = 20)
cat("done with model 8: binomial, corAR1, random effects, EddyDist, ts spline \n")


# ## try with random effects smooths
# cat("Try random effects smooths instead of thin plate \n")
# 
# gamm_full_TF$v09 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
#                         + s(SSH, bs="re", k=kVal)
#                         + s(log10_CHL, bs="re", k=kVal)
#                         + s(Neg_EddyDist, bs="re", k=kVal)
#                         + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                         + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                         + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                         + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                         data = transformedCovars.train,
#                         weights = myWeights,
#                         na.action = na.omit,family = nb(),
#                         control = list(opt='optim',maxIter = 100,msMaxIter=100),
#                         correlation = corAR1(form=~1|fac1),niterPQL = 20)
# cat("done with model 9: neg binom, corAR1, Neg_EddyDist, 're' spline \n")
# 
# gamm_full_TF$v10 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
#                         + s(SSH, bs="re", k=kVal)
#                         + s(log10_CHL, bs="re", k=kVal)
#                         + s(Neg_EddyDist, bs="re", k=kVal)
#                         + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                         + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                         + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                         + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                         data = transformedCovars.train,
#                         weights = myWeights,
#                         na.action = na.omit,family = binomial(),
#                         control = list(opt='optim',maxIter = 100, msMaxIter=100),
#                         correlation = corAR1(form=~1|fac1),niterPQL = 20)
# cat("done with model 10: binomial, corAR1, Neg_EddyDist, 're' spline \n")
# 
# gamm_full_TF$v11 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
#                         + s(SSH, bs="re", k=kVal)
#                         + s(log10_CHL, bs="re", k=kVal)
#                         + s(EddyDist, bs="re", k=kVal)
#                         + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                         + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                         + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                         + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                         data = transformedCovars.train,
#                         weights = myWeights,
#                         na.action = na.omit,family = nb(),
#                         control = list(opt='optim',maxIter = 100, msMaxIter=100),
#                         correlation = corAR1(form=~1|fac1),niterPQL = 20)
# cat("done with model 11: neg binom, corAR1, EddyDist, 're' spline \n")
# 
# gamm_full_TF$v12 <- gamm(y_TF~ s(SST, bs="re", k=kVal)
#                         + s(SSH, bs="re", k=kVal)
#                         + s(log10_CHL, bs="re", k=kVal)
#                         + s(EddyDist, bs="re", k=kVal)
#                         + s(log10_HYCOM_MAG_0,bs="re", k=kVal)
#                         + s(HYCOM_SALIN_0, bs="re", k=kVal)
#                         + s(log10_FrontDist_Cayula,bs="re", k=kVal)
#                         + s(log10_HYCOM_MLD, bs="re", k=kVal),
#                         data = transformedCovars.train,
#                         weights = myWeights,
#                         na.action = na.omit,family = binomial(),
#                         control = list(opt='optim',maxIter = 100,msMaxIter=100),
#                         correlation = corAR1(form=~1|fac1),niterPQL = 20)
# cat("done with model 12: binomial, corAR1, EddyDist, 're' spline \n")
# 

gamm_full_TF$v13 <- gamm(y_TF~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                weights = myWeights,
                                data = transformedCovars.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim',maxIter = 100,msMaxIter=100,niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 13: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_TF$v14 <- gamm(y_TF~ s(SST, bs="ts", k=kVal)
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                weights = myWeights,
                                data = transformedCovars.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim',maxIter = 100,msMaxIter=100, niterEM=0),
                                correlation = corAR1(form=~1|fac1))#
cat("done with model 14: quasibinomial, corAR1, EddyDist, ts spline \n")

gamm_full_TF$v15 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                weights = myWeights,
                                data = transformedCovars.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim',maxIter = 100,msMaxIter=100, niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 15: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")

gamm_full_TF$v16 <- gamm(y_TF~ s(SST, bs="ts", k=kVal) 
                                + s(SSH, bs="ts", k=kVal)
                                + s(log10_CHL, bs="ts", k=kVal)
                                + s(Neg_EddyDist, bs="ts", k=kVal)
                                + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                                + s(log10_FrontDist_Cayula,bs="ts", k=kVal)
                                + s(HYCOM_SALIN_0,bs="ts", k=kVal)
                                + s(log10_HYCOM_MLD, bs="ts", k=kVal),
                                weights = myWeights,
                                data = transformedCovars.train,
                                na.action = na.omit,family = quasibinomial(),
                                control = list(opt='optim',maxIter = 100,msMaxIter=100, niterEM=0),
                                correlation = corAR1(form=~1|fac1), random=list(fac1=~1))#
cat("done with model 16: quasibinomial, corAR1, random effects, EddyDist, ts spline \n")


model_AIC <- NULL
for (iMod in 1:length(gamm_full_TF)){
  model_AIC[iMod] = AIC(gamm_full_TF[[iMod]]$lme)
  cat('Model ', names(gamm_full_TF[iMod]), ' AIC = ', model_AIC[iMod], '\n')
}

best_Combined_model <- names(gamm_full_TF[which.min(model_AIC)])

cat("\n BEST COMBINED MODEL IS # ", best_Combined_model," \n")
cat("Saving all models to file \n \n")

save(gamm_full_TF, best_Combined_model, transformedCovars.train,
     transformedCovars.test,mergedTest.set,mergedTrain.set,
     file = paste(outDir,SP,'_AcousticAndVisual_binomial_GAMMs_ALL.Rdata',sep=''))



#################### Save best full model ################################################

gamm_pruned_TF_best <- gamm(y_TF~ s(SST, bs="ts", k=kVal)
                        + s(SSH, bs="ts", k=kVal)
                        + s(Neg_EddyDist, bs="ts", k=kVal)
                        + s(log10_HYCOM_MAG_0,bs="ts", k=kVal)
                        + s(HYCOM_SALIN_0, bs="ts", k=kVal)
                        + s(log10_FrontDist_Cayula,bs="ts", k=kVal),
                        data = transformedCovars.train,
                        weights = myWeights,
                        na.action = na.omit,family = nb(),
                        control = list(opt='optim',maxIter = 100,msMaxIter=100, niterEM=0,keepData = TRUE),
                        correlation = corAR1(form=~1|fac1),niterPQL = 20)

model <- gamm_pruned_TF_best$gam
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',
SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],
UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", 
                                  transformedCovars.train[,c("SST","SSH","Neg_EddyDist",
                                  "log10_HYCOM_MAG_0","log10_FrontDist_Cayula","HYCOM_SALIN_0")],
                                  NULL, y_TF, NULL, NULL, coordinateSystem, model)
# Then save 'model' to file (.Rdata)
save(model, modelMetadata, gamm_pruned_TF_best,
     file = paste(outDir,SP,'_gamm_pruned_TF_best.Rdata',sep='')) 

# Output best acoustic only model summary text to file
sink(paste0(outDir,SP,'_GAMM_TF_pruned_best_summary.txt'))
summary(model)
sink()

# Plot model smooths
png(paste0(outDir,SP,'_BestModel_smooths.png',sep=''),width = 2000,height = 1600,res=300)
plot(model,pages=1,cex.lab = 1.3,cex.axis = 1.1,scale=0)
dev.off()

# Predict on test data at MC, for comparison...
compSet <- which((mergedTest.set$fac1)==5)
pred <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars.test[compSet,],
                    type = 'response',na.action = na.omit)
png(paste0(outDir,SP,'_Combined_Prediction_timeseries.png',sep=''),width = 2000,height = 1500,res=300)
twoord.plot(mergedTest.set$date[compSet],pred,
            mergedTest.set$date[compSet],(mergedTest.set$Density[compSet]),xlab="Date",type=c("l","b"),
            main="Test Data vs. Model Prediction", ylab = "Model Probability of Occurrence",rylab = "Density From Data")
dev.off()

# Predict on acoustic test data, using combined modelfor comparison...
pred_MC <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MC,],
                       type = 'response',na.action = na.pass)
pred_GC <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_GC,],
                       type = 'response',na.action = na.pass)
pred_DT <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DT,],
                       type = 'response',na.action = na.pass)
pred_DC <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_DC,],
                       type = 'response',na.action = na.pass)
pred_MP <- predict.gam(gamm_pruned_TF_best$gam,transformedCovars_AcOnly.test[compAcSet_MP,],
                       type = 'response',na.action = na.pass)
occurIdx = which(as.POSIXct(pOccur[,1])>='2013-01-01' & as.POSIXct(pOccur[,1])<'2014-01-01')

png(paste0(outDir,SP,'_Full_Prediction_timeseries.png',sep=''),width = 2000,height = 3000,res=300)
par(mfrow = c(5,1))
twoord.plot(Test_AcOnly.set$date[compAcSet_MC],pred_MC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,2],type=c("l","l"),mar = c(1, 5, 1, 5),
            main="Test Data vs. Model Prediction",xtickpos = dateTicks,xticklab = rep("",5),axislab.cex=.8)
twoord.plot(Test_AcOnly.set$date[compAcSet_GC],pred_GC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,3],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DT],pred_DT,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,4],type=c("l","l"),
            ylab = "Model Probability of Occurrence",mar = c(1, 5, 0, 5),
            rylab = "Weekly Probability of Occurrence from Data",xtickpos = dateTicks,axislab.cex=.8,
            xticklab = rep("",5))
twoord.plot(Test_AcOnly.set$date[compAcSet_DC],pred_DC,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,5],type=c("l","l"),axislab.cex=.8,
            xtickpos =dateTicks,xticklab = rep("",5),mar = c(1, 5, 0, 5))
twoord.plot(Test_AcOnly.set$date[compAcSet_MP],pred_MP,
            as.POSIXct(pOccur[occurIdx,1]),pOccur[occurIdx,6],type=c("l","l"),mar = c(3, 5, 0, 5),
            xlab = "Date",xtickpos = dateTicks, xticklab = dateLabels,axislab.cex=.8)

#axis(1,at=dateTicks, labels=dateLabels)
dev.off()


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

