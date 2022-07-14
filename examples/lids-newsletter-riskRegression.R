### lids-newsletter-riskRegression.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Johan Sebastian Ohlendorff
## This R file illustrates in full details the R code shown in the newsletter. 
library(survival) # for using the pbc data set 
library(data.table)
library(riskRegression) # for using the Score function, synthesize function and various plot functions
library(lava) # necessary for the synthesize function
library(rms) # for using the rcs function
library(prodlim)
library(randomForestSRC) # for using random survival forests.
data(pbc,package="survival")
setDT(pbc)
pbc <- na.omit(pbc)
pbc$logbili <- log(pbc$bili)
set.seed(7)

# simulate training data based on real data (i.e. synthesize new data)
pbc_alike <- synthesize(object=Surv(time,status)~ sex + age + logbili + chol + hepato + spiders + protime + albumin + platelet + trig + trt + ast,data=pbc)
d <- simsynth(pbc_alike,n=400,seed=7)

# train the models on training data (uses splines in m1)
m1 <- CSC(Hist(time,status)~sex+rcs(age)+rcs(logbili)+rcs(protime)+hepato+rcs(chol)+spiders,data=d,fitter="cph")
m2 <- rfsrc(Surv(time,status)~sex+age+logbili+chol+hepato+spiders+protime+albumin+platelet+trig+trt+ast,data=d,fitter="cph",cause=2,seed=9)

# simulate test data and evaluate risk prediction models on test data 
nd <- simsynth(pbc_alike,3000,seed=28)
x <- Score(list("CSC"=m1,"FOREST"=m2),
           data=nd, formula=Hist(time,status)~1,
           cause=2, time=(0:5)*365.25,
           summary="risks",plots=c("calibration","Roc"))
summary(x, times=3*365.25, what="score")

# figures shown in newsletter 
plotRisk(x,times=3*365.25,models=c("FOREST","CSC"))
plotBrier(x)
plotAUC(x) # not shown
plotCalibration(x,times=3*365.25)
plotROC(x,times=3*365.25)

par(mfrow=c(2,1),mar=c(4,4,2,2))
setwd("~/research/Methods/LiDS/worg/")
plotRisk(x,times=3*365.25,models=c("FOREST","CSC"),legend.x=.58,legend.y=.5,legend.cex=.8,legend.title="Horizon: 3 years")
text("A",x=-.1,y=1.1,xpd=NA,cex=1.5)
# plotBrier(x,conf.int=FALSE,axis1.at=(0:5)*365.25,axis1.lab=0:5,xlab="Horizon (years)")
# text("B",x=-20,y=.33,xpd=NA,cex=1.5)
plotCalibration(x,times=3*365.25,models=c("FOREST","CSC"),auc.in.legend=0L,legend.title="Horizon: 3 years")
text("B",x=-.1,y=1.1,xpd=NA,cex=1.5)
# plotROC(x,times=3*365.25,models=c("FOREST","CSC"),legend.title="Horizon: 3 years")
# text("D",x=-.1,y=1.1,xpd=NA,cex=1.5)

par(mfrow=c(2,2),mar=c(4,4,2,2))
setwd("~/research/Methods/LiDS/worg/")
plotRisk(x,times=3*365.25,models=c("FOREST","CSC"),legend.x=.6,legend.y=.5,legend.title="Horizon: 3 years")
text("A",x=-.1,y=1.1,xpd=NA,cex=1.5)
plotBrier(x,conf.int=FALSE,axis1.at=(0:5)*365.25,axis1.lab=0:5,xlab="Horizon (years)")
text("B",x=-20,y=.33,xpd=NA,cex=1.5)
plotCalibration(x,times=3*365.25,models=c("FOREST","CSC"),auc.in.legend=0L,legend.title="Horizon: 3 years")
text("C",x=-.1,y=1.1,xpd=NA,cex=1.5)
plotROC(x,times=3*365.25,models=c("FOREST","CSC"),legend.title="Horizon: 3 years")
text("D",x=-.1,y=1.1,xpd=NA,cex=1.5)

