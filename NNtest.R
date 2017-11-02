# Set a seed
set.seed(600)

# make fac4 for training
fac4 <- array(data = NA, dim = c(length(transformedCovars.train$fac1),1))
weightsG0<- array(data = 1, dim = c(length(transformedCovars.train$fac1),1))
for (iFac in 1:length(transformedCovars.train$fac1)) {
  if (!is.na(transformedCovars.train$fac1[iFac])){
    if (transformedCovars.train$fac1[iFac]>5) {
      fac4[iFac,1] <- 0
      if (mergedTrain.set$Density[iFac]==0){
        # if it's visual data and it's a zero, adjust by g0 ie, only a 30% chance it was a true zero.
        weightsG0[iFac,1] <- .3
        
      }
    } else {
      fac4[iFac,1] <- transformedCovars.train$fac1[iFac]
    }
  }
}
transformedCovars.train$fac4<-fac4
transformedCovars.train$weightsG0<-weightsG0

# make fac4 for testing
fac4 <- array(data = NA, dim = c(length(transformedCovars.test$fac1),1))
weightsG0<- array(data = 1, dim = c(length(transformedCovars.test$fac1),1))
for (iFac in 1:length(transformedCovars.test$fac1)) {
  if (!is.na(transformedCovars.test$fac1[iFac])){
    if (transformedCovars.test$fac1[iFac]>5) {
      # it's visual data
      fac4[iFac,1] <- 0
      if (mergedTrain.set$Density[iFac]==0){
        # if it's visual data and it's a zero, adjust by g0 ie, only a 30% chance it was a true zero.
        weightsG0[iFac,1] <- .3
        }
    } else {
      fac4[iFac,1] <- transformedCovars.test$fac1[iFac]
      weightsG0[iFac,1] <- 1
    }
  }
}
transformedCovars.test$fac4<-fac4
transformedCovars.test$weightsG0<-weightsG0

# Set up training data
allDataTrain <- transformedCovars.train
allDataTrain$y <- mergedTrain.set$Density
allDataTrain$y_TF <- (mergedTrain.set$Density>0)*1

goodData <- which(!is.na(rowSums(allDataTrain)))
trainDataNoNa <- allDataTrain[goodData,]
apply(trainDataNoNa,2,function(x) sum(is.na(x)))


allDataTest <- transformedCovars.test
allDataTest$y <- mergedTest.set$Density
allDataTest$y_TF <-  (mergedTest.set$Density>0)*1

goodData <- which(!is.na(rowSums(allDataTest)))
testDataNoNa <- allDataTest[goodData,]

# Remove missing data
apply(testDataNoNa,2,function(x) sum(is.na(x)))

library(neuralnet)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Neural net fitting

# Scaling training data for the NN
maxsTrain <- apply(trainDataNoNa, 2, max) 
minsTrain <- apply(trainDataNoNa, 2, min)
train_scaledData <- as.data.frame(scale(trainDataNoNa, center = minsTrain, scale = maxsTrain-minsTrain))
# Scaling test data for the NN
maxsTest <- apply(testDataNoNa, 2, max) 
minsTest <- apply(testDataNoNa, 2, min)
test_scaledData <- as.data.frame(scale(testDataNoNa, center = minsTest, scale = maxsTest-minsTest))

# NN training
n <- names(train_scaledData)

sigmoid = function(x) {
  1 / (1 + exp(-x))
}
myIndices <- c(1:5,7,10)
nMax <- length(myIndices)

f <- as.formula(paste("y_TF ~", paste(n[myIndices], collapse = " + ")))
nn <- neuralnet(f,data=train_scaledData,hidden=c(6,2), lifesign = "full",threshold = 0.03,stepmax = 1e+05,
                linear.output=FALSE, act.fct = "logistic", rep = 1)
nn2 <- nnet(f, data=train_scaledData, size=5, na.action = na.omit,  rang = 0.1,weights= trainDataNoNa$weightsG0,
            decay = 2e-4, maxit = 1000)


# Predict
pr.nn <- compute(nn,test_scaledData[,myIndices])
pr.nn.train <- compute(nn,train_scaledData[,myIndices])
pr.nn2 <- predict(nn2,test_scaledData[,myIndices])
pr.nn2.train <- predict(nn2,train_scaledData[,myIndices])


# Results from NN are normalized (scaled)
# Descaling for comparison
pr.nn_ <- pr.nn$net.result*(max(trainDataNoNa$y_TF)-min(trainDataNoNa$y_TF))+min(trainDataNoNa$y_TF)
test.r <- (testDataNoNa$y_TF)*(max(testDataNoNa$y_TF)-min(testDataNoNa$y_TF))+min(testDataNoNa$y_TF)

pr.nn2_ <- pr.nn2*(max(trainDataNoNa$y_TF)-min(trainDataNoNa$y_TF))+min(trainDataNoNa$y_TF)


# Descaling for comparison
pr.nn.train_ <- pr.nn.train$net.result*(max(trainDataNoNa$y_TF)-min(trainDataNoNa$y_TF))+min(trainDataNoNa$y_TF)
train.r <- (trainDataNoNa$y_TF)*(max(trainDataNoNa$y_TF)-min(trainDataNoNa$y_TF))+min(trainDataNoNa$y_TF)

pr.nn2.train_ <- pr.nn2.train*(max(trainDataNoNa$y_TF)-min(trainDataNoNa$y_TF))+min(trainDataNoNa$y_TF)


# Calculating MSE
MSE.nn.test <- sum((test.r - pr.nn_)^2)/nrow(test_scaledData[1:nMax])
MSE.nn.train <- sum((train.r - pr.nn.train_)^2)/nrow(train_scaledData[1:nMax])
MSE.nn.test
MSE.nn.train

MSE.nn2.test <- sum((test.r - pr.nn2_)^2)/nrow(test_scaledData[1:nMax])
MSE.nn2.train <- sum((train.r - pr.nn2.train_)^2)/nrow(train_scaledData[1:nMax])
MSE.nn2.test
MSE.nn2.train

# Visual plot of the model
plot(nn)

plot(test.r,pr.nn_ ,col='red',main='Real vs predicted NN: TEST SET',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(train.r,pr.nn.train_ ,col='red',main='Real vs predicted NN: TRAIN SET',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')