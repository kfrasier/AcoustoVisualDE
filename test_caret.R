#AcOnly_train_scaled$yAcOnly_TF <- make.names(factor(AcOnly_train_scaled$yAcOnly_TF))
f.AcOnly_NN1 <- as.formula(paste("yAcOnly_TF ~", paste(n[model1.indices], collapse = " + ")))

fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3)# classProbs = TRUE

nnetGrid <-  expand.grid(size = seq(from = 1, to = 15, by = 2),
                         decay = [0.01, .001,.0001],bag = TRUE)


t1<-train(f.AcOnly_NN1, data=AcOnly_train_scaled,method = 'avNNet',
      metric="RMSE",na.action = na.omit,
      trControl = fitControl,
      tuneGrid = nnetGrid)



t2<-train(f.AcOnly_NN1, data=AcOnly_train_scaled,method = 'avNNet',
          metric="RMSE",na.action = na.omit,
          trControl = fitControl,
          tuneGrid = nnetGrid,importance = TRUE)