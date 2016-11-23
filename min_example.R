min.example<-function(ddfData,keyList){
  
  for (i1 in 1:length(keyList)){

    detFun1[[i1]] <-ddf(method ='ds', dsmodel =~ mcds(key = as.character(keyList[i1]), formula = ~ 1),
                        data = as.data.frame(ddfData), meta.data = list(binned=F, width=5000, left=0))
    
  }
  
}