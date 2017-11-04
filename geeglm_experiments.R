library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(yags)            # for the GEEs (QIC provided)
library(splines)         # to construct the B-splines within a GEE-GLM
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots

allDataTrain <- transformedCovars_VisOnly.train
allDataTrain$yVisOnly <-  yVisOnly
allDataTrain$yVisOnlyTF <-  yVisOnly_TF

goodData <- which(!is.na(rowSums(allDataTrain)))
dataNoNa <- allDataTrain[goodData,]


yVisOnlyTF <- dataNoNa$yVisOnlyTF

modelFull1 <- geeglm(yVisOnlyTF ~ bs(SST, knots = mean(SST)) + bs(SSH, knots = mean(SSH))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula))+
                     bs(EddyDist, knots = mean(EddyDist)),
                     family = binomial, data = dataNoNa, corstr = "ar1",id = fac2)


modelFull1_z <- zeroinfl(yVisOnlyTF ~ bs(SST, knots = mean(SST))+bs(SSH, knots = mean(SSH))+bs(log10_CHL, knots = mean(log10_CHL)),
                         data = dataNoNa, dist = 'poisson')

q1 <- qic(modelFull1)$QIC
modelFull2 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + bs(SSH, knots = mean(SSH))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula))+
                     bs(Neg_EddyDist, knots = mean(Neg_EddyDist)),
                     family = binomial, data = dataNoNa, corstr = "ar1",id = fac2) 
q2 <- qic(modelFull2)$QIC
modelFull3 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + bs(SSH, knots = mean(SSH))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula))+
                     bs(Pos_EddyDist, knots = mean(Pos_EddyDist)),
                     family = binomial, data = dataNoNa, corstr = "ar1",id = fac2) 
q3 <- qic(modelFull3)$QIC
modelFull4 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                     family = binomial, data = dataNoNa, corstr = "ar1",id = fac2) 
q4 <- qic(modelFull4)$QIC

model<-c("full_Eddy","full_negEddy","full_posEddy","full_SSHandEddy")
QIC<-c(q1,q2,q3,q4)
data.frame(rbind(model,QIC))

# Lowest QIC model is interaction btwn SSH and Eddy

# Try linear terms instead of splines 
# SST is better as linear
modelFull4_1 <- geeglm(yAcOnlyTF ~ SST + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                    family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)  # 3260.536601
# Linear SSH, no improvement
modelFull4_2 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (SSH*bs(EddyDist, knots = mean(EddyDist)))+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                     family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) 
# Linear eddyDist, no improvement
modelFull4_3 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*EddyDist)+
                     bs(log10_CHL, knots = mean(log10_CHL)) + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                     family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) 
# Linear log10_CHL, no improvement
modelFull4_4 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                     log10_CHL + bs(log10_HYCOM_MLD, knots = mean(log10_HYCOM_MLD))+
                     bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                     bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                     family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) 
# Linear MLD is better, keep it
modelFull4_5 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                       bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                      family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) # 3255.621459
# Linear Salinity, no improvement
modelFull4_6 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       HYCOM_SALIN_0 + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                       bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                       family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) 
# Linear Mag, no improvement
modelFull4_7 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + log10_HYCOM_MAG_0 +
                       bs(HYCOM_UPVEL_50, knots = mean(HYCOM_UPVEL_50)) + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                       family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) 
# Linear HYCOM_UPVEL_50 is better, keep it
modelFull4_8 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                       HYCOM_UPVEL_50 + bs(log10_FrontDist_Cayula, knots = mean(log10_FrontDist_Cayula)),
                       family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) # 3251.101557
# Linear log10_FrontDist_Cayula is better, keep it
modelFull4_9 <- geeglm(yAcOnlyTF ~ bs(SST, knots = mean(SST)) + (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                       HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                       family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) # 3246.719926

# Now start deleting variables:
# pull SST out, slightly improved
modelReduced_1 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                       bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                       bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                       HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                       family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)# 3244.776771
# remove SSH -> worse
modelReduced_2 <- yags(yAcOnlyTF ~ bs(EddyDist, knots = mean(EddyDist))+
                         bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)
# remove EddyDist -> worse
modelReduced_3 <- yags(yAcOnlyTF ~ bs(SSH, knots = mean(SSH))+
                         bs(log10_CHL, knots = mean(log10_CHL)) + log10_HYCOM_MLD +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)
# remove log10_CHL -> worse
modelReduced_4 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         log10_HYCOM_MLD +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)
# remove log10_HYCOM_MLD -> better
modelReduced_5 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         bs(log10_CHL, knots = mean(log10_CHL)) +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0)) + bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)# 3243.870873
# remove HYCOM_SALIN_0 -> worse
modelReduced_6 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         bs(log10_CHL, knots = mean(log10_CHL)) +
                         bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)
# remove log10_HYCOM_MAG_0 -> worse
modelReduced_7 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         bs(log10_CHL, knots = mean(log10_CHL)) +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0))+
                         HYCOM_UPVEL_50 + log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1)
# remove HYCOM_UPVEL_50 -> better
modelReduced_8 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         bs(log10_CHL, knots = mean(log10_CHL)) +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0))+ bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0))+
                         log10_FrontDist_Cayula,
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) #3241.148521
# remove log10_FrontDist_Cayula -> better
modelReduced_9 <- yags(yAcOnlyTF ~ (bs(SSH, knots = mean(SSH))*bs(EddyDist, knots = mean(EddyDist)))+
                         bs(log10_CHL, knots = mean(log10_CHL)) +
                         bs(HYCOM_SALIN_0, knots = mean(HYCOM_SALIN_0))+ bs(log10_HYCOM_MAG_0, knots = mean(log10_HYCOM_MAG_0)),
                         family = binomial, data = dataNoNa, corstruct = "independence",id = fac1) #3237.444515
