library(survival)
library(timeROC)
##-------------Without competing risks-------------------
data(pbc)
head(pbc)
pbc<-pbc[!is.na(pbc$trt),] # select only randomised subjects
pbc$status<-as.numeric(pbc$status==2) # create event indicator: 1 for death, 0 for censored     
# Se, Sp, PPV and NPV computation for serum bilirunbin at threshold c=0.9(mg/dl) 
res.SeSpPPVNPV.bili <- SeSpPPVNPV(cutpoint=0.9,
                                  T=pbc$time,
                                  delta=pbc$status,marker=pbc$bili,
                                  cause=1,weighting="marginal",
                                  times=1200,
                                  iid=TRUE)
x<-Score(list(pbc$bili),formula=Surv(time,status)~1,data=pbc,times=1200,metrics="auc", cutpoints=0.9)
## Check that they give the same results
x$AUC$res.cut[,c(4,6,8,10)]
res.SeSpPPVNPV.bili[c(1,2,3,4)]
x$AUC$res.cut[,c(5,7,9,11)]
res.SeSpPPVNPV.bili$inference[c(7,9,10,12)]

## Binary
res.SeSpPPVNPV.bili <- SeSpPPVNPV(cutpoint=0.9,
                                  T=pbc$time,
                                  delta=rep(1,nrow(pbc)),marker=pbc$bili,
                                  cause=1,weighting="marginal",
                                  times=1200,
                                  iid=TRUE)
pbc$Y <- 1*(pbc$time <= 1200)
x<-Score(list(pbc$bili),formula=Y~1,data=pbc,metrics="auc", cutpoints=0.9)
## Check that they give the same results
x$AUC$res.cut[,c(2,4,6,8)]
res.SeSpPPVNPV.bili[c(1,2,3,4)]
x$AUC$res.cut[,c(3,5,7,9)]
res.SeSpPPVNPV.bili$inference[c(7,9,10,12)]

##-------------With competing risks-------------------
## Also test computation time
data(Paquid)
system.time(res.SeSpPPVNPV.DSST <- SeSpPPVNPV(cutpoint=22,
                                  T=Paquid$time,
                                  delta=Paquid$status,marker=Paquid$DSST,
                                  cause=1,weighting="marginal",
                                  times=5,iid=TRUE))

library(riskRegression)
system.time(x<-Score(list(Paquid$DSST),formula=Hist(time,status)~1,data=Paquid,times=c(5),metrics="auc", cutpoints=22))
## Check that they give the same results
x$AUC$res.cut[,c(4,6,8,10)]
res.SeSpPPVNPV.DSST[c(1,3,4,6)]
x$AUC$res.cut[,c(5,7,9,11)]
res.SeSpPPVNPV.DSST$inference[c(7,9,10,12)]

## Test Cox censoring
x<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
x$AUC$res.cut
y<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22,censoring.save.memory = TRUE)
y$AUC$res.cut
## They are a bit different.
