# quantiles
qASST<-quantile(c(Train_AcOnly.set2$SST,Test_AcOnly.set2$SST),c(.05,.95))
qVSST<-quantile(c(Train_VisOnly.set2$SST,Test_VisOnly.set2$SST),c(.05,.95))

qASSH<-quantile(c(Train_AcOnly.set2$SSH,Test_AcOnly.set2$SSH),c(.05,.95))
qVSSH<-quantile(c(Train_VisOnly.set2$SSH,Test_VisOnly.set2$SSH),c(.05,.95))

qACHL<-quantile(c(Train_AcOnly.set2$CHL,Test_AcOnly.set2$CHL),c(.05,.95), na.rm=TRUE)
qVCHL<-quantile(c(Train_VisOnly.set2$CHL,Test_VisOnly.set2$CHL),c(.05,.95), na.rm=TRUE)

qASAL<-quantile(c(Train_AcOnly.set2$HYCOM_SALIN_0,Test_AcOnly.set2$HYCOM_SALIN_0),c(.05,.95))
qVSAL<-quantile(c(Train_VisOnly.set2$HYCOM_SALIN_0,Test_VisOnly.set2$HYCOM_SALIN_0),c(.05,.95))

qAMLD<-quantile(c(Train_AcOnly.set2$HYCOM_MLD,Test_AcOnly.set2$HYCOM_MLD),c(.05,.95))
qVMLD<-quantile(c(Train_VisOnly.set2$HYCOM_MLD,Test_VisOnly.set2$HYCOM_MLD),c(.05,.95))

qAUPVEL<-quantile(c(Train_AcOnly.set2$HYCOM_UPVEL_50,Test_AcOnly.set2$HYCOM_UPVEL_50),c(.05,.95))
qVUPVEL<-quantile(c(Train_VisOnly.set2$HYCOM_UPVEL_50,Test_VisOnly.set2$HYCOM_UPVEL_50),c(.05,.95))

qAEddy<-quantile(c(Train_AcOnly.set2$EddyDist,Test_AcOnly.set2$EddyDist),c(.05,.95))
qVEddy<-quantile(c(Train_VisOnly.set2$EddyDist,Test_VisOnly.set2$EddyDist),c(.05,.95))

qAMAG<-quantile(c(Train_AcOnly.set2$HYCOM_MAG_0,Test_AcOnly.set2$HYCOM_MAG_0),c(.05,.95))
qVMAG<-quantile(c(Train_VisOnly.set2$HYCOM_MAG_0,Test_VisOnly.set2$HYCOM_MAG_0),c(.05,.95))

# means
mASST<-mean(c(Train_AcOnly.set2$SST,Test_AcOnly.set2$SST))
mVSST<-mean(c(Train_VisOnly.set2$SST,Test_VisOnly.set2$SST))

mASSH<-mean(c(Train_AcOnly.set2$SSH,Test_AcOnly.set2$SSH))
mVSSH<-mean(c(Train_VisOnly.set2$SSH,Test_VisOnly.set2$SSH))

mACHL<-mean(c(Train_AcOnly.set2$CHL,Test_AcOnly.set2$CHL), na.rm=TRUE)
mVCHL<-mean(c(Train_VisOnly.set2$CHL,Test_VisOnly.set2$CHL), na.rm=TRUE)

mASAL<-mean(c(Train_AcOnly.set2$HYCOM_SALIN_0,Test_AcOnly.set2$HYCOM_SALIN_0))
mVSAL<-mean(c(Train_VisOnly.set2$HYCOM_SALIN_0,Test_VisOnly.set2$HYCOM_SALIN_0))

mAMLD<-mean(c(Train_AcOnly.set2$HYCOM_MLD,Test_AcOnly.set2$HYCOM_MLD))
mVMLD<-mean(c(Train_VisOnly.set2$HYCOM_MLD,Test_VisOnly.set2$HYCOM_MLD))

mAUPVEL<-mean(c(Train_AcOnly.set2$HYCOM_UPVEL_50,Test_AcOnly.set2$HYCOM_UPVEL_50))
mVUPVEL<-mean(c(Train_VisOnly.set2$HYCOM_UPVEL_50,Test_VisOnly.set2$HYCOM_UPVEL_50))

mAEddy<-mean(c(Train_AcOnly.set2$EddyDist,Test_AcOnly.set2$EddyDist))
mVEddy<-mean(c(Train_VisOnly.set2$EddyDist,Test_VisOnly.set2$EddyDist))

mAMAG<-mean(c(Train_AcOnly.set2$HYCOM_MAG_0,Test_AcOnly.set2$HYCOM_MAG_0))
mVMAG<-mean(c(Train_VisOnly.set2$HYCOM_MAG_0,Test_VisOnly.set2$HYCOM_MAG_0))

table(cbind(rbind(mASST,mVSST,mASSH,mVSSH,mACHL,mVCHL,
            mASAL,mVSAL,mAMLD,mVMLD,mAUPVEL,mVUPVEL,mAEddy,
            mVEddy,mAMAG,mVMAG),
  rbind(qASST,qVSST,qASSH,qVSSH,qACHL,qVCHL,
            qASAL,qVSAL,qAMLD,qVMLD,qAUPVEL,qVUPVEL,qAEddy,
            qVEddy,qAMAG,qVMAG)))
        

