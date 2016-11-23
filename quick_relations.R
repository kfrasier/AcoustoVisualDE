t2 <-covarListTrans.train[,c(1,4,10)]


plot.covarDensity.proposal(t2,colnames(t2), mod.train.set$SpPresent)

c("SST","log_10(Mixed Layer Depth)", "Dist. to Eddy")


encounterAll <- gam(y ~ s(SST_daily, bs="ts",k=5)+ s(SSH_daily,bs="ts",k=5) 
                    + s(log10_CHL_daily_climate,bs="ts",k=5) ,
                    method = "GCV.Cp", data = transformedCovars.train, family = tw(),
                    offset = log(mergedTrain.set$EffectiveArea),na.action = na.omit)# 