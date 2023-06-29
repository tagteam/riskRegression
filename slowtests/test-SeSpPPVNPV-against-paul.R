##-------------Without competing risks-------------------
library(survival)
library(timeROC)
data(Paquid)
system.time(res.SeSpPPVNPV.DSST <- SeSpPPVNPV(cutpoint=22,
                                  T=Paquid$time,
                                  delta=Paquid$status,marker=Paquid$DSST,
                                  cause=1,weighting="marginal",
                                  times=5,iid=TRUE))

library(riskRegression)
system.time(x<-Score(list(Paquid$DSST),formula=Hist(time,status)~1,data=Paquid,times=c(5),metrics="auc", cutpoints=22))
x$AUC$res.cut
res.SeSpPPVNPV.DSST[c(1,3,4,6)]
res.SeSpPPVNPV.DSST$inference[c(7,9,10,12)]

## Test Cox censoring
x<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
x$AUC$res.cut
y<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22,censoring.save.memory = TRUE)
y$AUC$res.cut
## They are a bit different.
