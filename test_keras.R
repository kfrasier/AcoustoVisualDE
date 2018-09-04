
library(keras)

# AcOnly_train_scaled[,model1.indices]
#Multiple linear regression
modelAc <- keras_model_sequential() 
modelAc %>% 
  layer_dense(units = 64, input_shape = 8, activation='linear') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 64, input_shape = 8, activation='linear') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1, input_shape = 64, activation='linear') %>%
  layer_dropout(rate = 0.5) %>%
  summary(modelAc)

modelAc %>% compile(
  loss = 'mse',
  optimizer = 'adam',
  metrics = c("mean_absolute_percentage_error")
)

myY <- AcOnly_train_scaled$yAcOnlySqrt
yNoNA <- which(!is.na(myY))

myY<- as.numeric(AcOnly_train_scaled$yAcOnlySqrt[yNoNA])
myData <- as.matrix(AcOnly_train_scaled[yNoNA,model1.indices])


historylr <- modelAc %>% fit(
  myData, myY, 
  epochs = 20, batch_size = 100,shuffle = TRUE
)

myTestY<- AcOnly_test_scaled$yAcOnlySqrt
TestyNoNA <- which(!is.na(myTestY))
myTestY <- AcOnly_test_scaled$yAcOnlySqrt[TestyNoNA]
myTestData <-  as.matrix(AcOnly_test_scaled[TestyNoNA,model1.indices])
results <- modelAc %>% evaluate(myTestData, myTestY)


trainPredict<-predict_proba(modelAc,myData)
trainPredict[trainPredict<0]<-0
plot(trainPredict)
MSEtrain <- sum((trainPredict-AcOnly_train_scaled$yAcOnlySqrt)^2)/length(trainPredict)


testPredict<-predict_proba(modelAc,myTestData)
testPredict[testPredict<0]<-0
plot(testPredict)
MSEtest <- sum((testPredict-AcOnly_test_scaled$yAcOnlySqrt)^2)/length(testPredict)
