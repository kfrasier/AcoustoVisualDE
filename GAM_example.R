# Acoustovisual density estimation
library(mrds)
library(psych)
library(lubridate)
library(magic)
library(mgcv)

## modify these to match your paths
setwd("E:\\NASData")
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_missingdata.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_cleveland.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/transform_covars.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/plot_covarDensity.R')
source('E:/NASData/AcoustoVisualDE/AcoustoVisualDE/GetModelMetadata.R')


#### Parameters needed:  ####
outDir <- "E:/NASData/ModelData/" 
## Detection files
acousticSegFile <- "E:/NASData/ALLSITES_Ssp_density8hrs_jahStart.csv" # acoustic input file


############### Begin Calculations ################
graphics.off()
closeAllConnections()

## Load & Prune Acoustic data
cat("Loading acoustic data \n")


acSegmentsCSV <- read.csv(acousticSegFile, header = TRUE,na.strings=c(""," ","NA","-99999","NaN"))
acSegments <- as.data.frame(acSegmentsCSV)

## PRobably need to help it interpret the date field correctly:
acSegments$date <- c(as.Date(acSegments,"%m/%d/%Y")) #date

# Get the column names of the predictor variables. Modify the indices to select the ones you want
covarList<-names(acSegments[2:length(acSegments)])

###################################
# Oceanographic variables

# Explore the data, graphical output
cat("Begin exploratory plot generation\n")

# histograms of missing data
percFilled <- plot.missingdata(acSegments,covarList) 

# If you decide from the missing data plots that you want to restrict years going forward:
yearListIdx = as.numeric(format(acSegments$date,"%Y"))

# separate into training and test sets:
keepDates.train <- which(yearListIdx != 2009 & yearListIdx >= 2003 & yearListIdx <= 2012)
keepDates.test <- which(yearListIdx == 2009)

acSegmentsTrainSet<- acSegments[keepDates.train,]
acSegmentsTestSet<- acSegments[keepDates.test,]
  

# Cleveland dot plots:
# Look at the data with no tranforms to see if it needs them
plot.cleveland(acSegmentsTrainSet,covarList,FALSE)


# Choose set to keep
covarList2 <- c("Density","SST","FrontDist_Cayula")
# and choose transforms
transformList <- c("none","none","log10")


# restrict covariates again to limited set
acSegmentsTrainSet2<- acSegmentsTrainSet[,covarList2]
acSegmentsTestSet2<- acSegmentsTestSet[,covarList2]


transformedCovars.train <- transform.covars(acSegmentsTrainSet2,covarList2,transformList)
transformedCovars.test <- transform.covars(acSegmentsTestSet2,covarList2,transformList)

# Now plot with transformations
plot.cleveland(transformedCovars.train,colnames(transformedCovars.train),TRUE)

# presence absence histograms - this only works with presence absence data
# plot.covarDensity(transformedCovars.train,colnames(transformedCovars.train),mergedTrain.set$SpPresent)

# correlation
png(paste(outDir,SP,'_correlations_noTransform.png',sep=''), width = 2000, height = 1600)
pairs.panels(mergedTrain.set2, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off()

png(paste(outDir,SP,'_correlations_withTransform.png',sep=''), width = 2000, height = 1600)
pairs.panels(transformedCovars.train, ellipses=FALSE, method = "spearman",cex.cor=.75)
dev.off() 

cat("Exploratory plots done\n")

###########################
# Run & evaluate models

y <- (acSegmentsTrainSet$MyResponseVariable)
plot(y)

cat("Run full GAM on with shrinkage\n")

# here's a GAM with all the possible predictor variables I had.
# I used a tweedie distribution ( family = tw() ), because my data are zero-inflated, but you may want a different one
# Note circular smooth (bs="cc") for polar variables like HYCOM_DIR (current direction, 0-360 deg)
encounterAll <- gam(y ~ s(SST, bs="ts",k=5) + s(SSH, bs="ts",k=5) + s(log10_CHL, bs="ts",k=5) + s(log10_HYCOM_MLD, bs="ts",k=5) + s(HYCOM_EMP, bs="ts",k=5)
         + s(HYCOM_DIR_0, bs="cc",k=5) + s(HYCOM_DIR_100, bs="cc",k=5) + s(HYCOM_DIR_800, bs="cc",k=5)
         + s(HYCOM_SALIN_0, bs="ts",k=5) + s(HYCOM_SALIN_100, bs="ts",k=5) + s(HYCOM_SALIN_800, bs="ts",k=5)
         + s(log10_HYCOM_MAG_0, bs="ts",k=5) + s(log10_HYCOM_MAG_100, bs="ts",k=5) + s(HYCOM_MAG_800, bs="ts",k=5)
         + s(HYCOM_UPVEL_100, bs="ts",k=5) + s(HYCOM_UPVEL_800),
         method = "REML", data = transformedCovars.train, family = tw(),
         na.action = na.omit)# 

# Output summary text to file
sink(paste(outDir,SP,'_GAM_density_full.txt'))
summary(encounterAll$gam)
gam.check(encounterAll$gam)
sink()

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll)
png(paste(outDir,SP,'_residuals_density_Full.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 

plot(encounterAll$gam,pages=1,ylim =c(-2,2))

cat("Run reduced GAM with without unused variables\n")
# look at that output, some covariates have been shrunk down to nothing, so remove them
ctl <- gam.control()
ctl$keepData = TRUE
encounterAll <- gam(y ~ s(SSH_daily, bs="ts",k=5)
                    + s(Month, bs="cc", k=5) 
                    + s(SST_daily, bs="ts", k=5)
                    + s(HYCOM_dir_daily, bs="cc",k=5) ,
                    control = list(keepData=TRUE),
                    method = "RML", data = transformedCovars.train, family = tw(),
                    na.action = na.omit)# 
# You could copy some of the plotting steps from above down to here to see the same plots for this.


# There should be a step after the modeling where you compare with the test set to evaluate prediction power.
# Haven't implemented that yet.



##############################################
# This last part is for making a prediction for the test year, or something, I've left it here, 
# but talk to me if you want to get it working because it's not trivial.

model <- encounterAll
coordinateSystem <- "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID[' GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
modelMetadata <- GetModelMetadata(terms(model), "mgcv", transformedCovars.train, NULL, y, NULL, NULL, coordinateSystem, model)
  
plot(encounterAll,pages=1,ylim =c(-2,2))

# Output summary text to file
sink(paste(outDir,SP,'_GAM_density_pruned.txt'))
summary(encounterAll$gam)
gam.check(encounterAll$gam)
sink()

plot(encounterAll$gam)

# Calculate and save residuals to text file
rsd <-residuals.gam(encounterAll$gam)
png(paste(outDir,SP,'_residuals_density_pruned.png',sep=''), width = 1000, height = 800)
plot(rsd)
dev.off() 

# save model
save(model, modelMetadata, file = paste(outDir,SP,'_GAM_density_pruned.Rdata',sep=''))

# Predict on test data
yTest <- mergedTest.set$SpEncounter_g0adj
pred <- predict.gam(encounterAll,transformedCovars.test, type = 'response',na.action = na.omit)
plot(yTest,pred)
