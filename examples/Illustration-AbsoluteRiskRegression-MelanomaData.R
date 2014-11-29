## load libraries and data

library(riskRegression)
library(pec)    ## to estimate the prediction performance
library(cmprsk) ## for comparison with the Fine-Gray model
data(Melanoma)
Melanoma$logthick <- log(Melanoma$thick)
help(Melanoma)

## censoring distribution

coxph(Surv(time,status==0)~sex+epicel+ulcer+age+logthick,data=Melanoma)

## Absolute risk regression
arr.invasion <- ARR(Hist(time,status)~invasion,
                    data=Melanoma,
                    cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
print(arr.invasion)

## Compare with stratified Aalen-Johansen
aj.invasion <- prodlim(Hist(time,status)~invasion,data=Melanoma)
plot(aj.invasion,cause=1,confint=FALSE)
plot(arr.invasion,newdata=data.frame(invasion=levels(Melanoma$invasion)),add=TRUE,lty=3,col=1:3)


## Fine-Gray regression model
library(cmprsk)
fg.invasion <- FGR(Hist(time,status)~invasion,data=Melanoma,cause=1)
print(fg.invasion)

## combined cause-specific Cox regression models (same formula for both causes)
## see below for an example with different formulae 
cox.invasion <- CSC(Hist(time,status)~invasion,data=Melanoma)
print(cox.invasion)

## Absolute risk regression with continuous predictor
arr.thick <- ARR(Hist(time,status)~thick,
                 data=Melanoma,
                 cause=1,
                 cens.model="cox",
                 cens.formula=~sex+epicel+ulcer+age+logthick)
plot(arr.thick)

## Nearest neighbor Aalen-Johansen at quantiles of thickness
aj.thick <- prodlim(Hist(time,status)~thick,data=Melanoma)
plot(aj.thick,cause=1,newdata=data.frame(thick=quantile(Melanoma$thick)),confint=FALSE,xlim=c(0,3000))
plot(arr.thick,newdata=data.frame(thick=quantile(Melanoma$thick)),add=TRUE,lty=3,col=1:5)

## putting thickness on log-scale
Melanoma$logthick <- log(Melanoma$thick)
arr.logthick <- ARR(Hist(time,status)~logthick,
                    data=Melanoma,
                    cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
plot(aj.thick,cause=1,newdata=data.frame(thick=quantile(Melanoma$thick)),confint=FALSE,xlim=c(0,3000))
plot(arr.logthick,newdata=data.frame(logthick=quantile(Melanoma$logthick)),add=TRUE,lty=3,col=1:5)

## comparison to logistic risk model
lrr.logthick <- LRR(Hist(time,status)~logthick,
                    data=Melanoma,
                    cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
plot(lrr.logthick)
plot(arr.logthick,newdata=data.frame(logthick=quantile(Melanoma$logthick)),add=TRUE,lty=3,col=1:5)

## comparison to combined cause-specific Cox regression models
csc.logthick <- CSC(Hist(time,status)~logthick,data=Melanoma,cause=1)
plot(csc.logthick)
plot(arr.logthick,newdata=data.frame(logthick=quantile(Melanoma$logthick)),add=TRUE,lty=3,col=1:5)

## Absolute risk regression
arr.fit <- ARR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
print(arr.fit)

## Logistic risk regression
lrr.fit <- LRR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
print(lrr.fit)

## Fine-Gray regression model
fg.fit <- FGR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1)
print(fg.fit)

## combined cause-specific Cox models
cox.fit <- CSC(list(cause1=Hist(time,status)~sex+epicel+ulcer+age+logthick,
                    cause2=Hist(time,status)~sex+age),data=Melanoma)
print(cox.fit)

# ARR versus combined Cox regression
# after 5 years = 1826.25 days
plot(predictEventProb(arr.fit,times=1826.25,newdata=Melanoma),
     predictEventProb(cox.fit,times=1826.25,newdata=Melanoma),
     main="Predicted probabilities of death due to malignant melanoma\nafter five years",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Absolute risk model",
     ylab="Combined cause-specific Cox models",
     axes=FALSE)
PercentAxis(1,at=seq(0,1,.25))
PercentAxis(2,at=seq(0,1,.25))
abline(a=0,b=1)


## prediction error 
perr.link <- pec(list(AbsRisk=arr.fit,CauseSpecCox=cox.fit,FineGray=fg.fit,LogisticRisk=lrr.fit),
                  formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,
                  cens.model="cox",
                  data=Melanoma,
                  maxtime=2500)
print(perr.link)

plot(perr.link)

cv.perr.link <- pec(list(AbsRisk=arr.fit,CauseSpecCox=cox.fit,FineGray=fg.fit,LogisticRisk=lrr.fit),
                  formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,
                  cens.model="cox",
                   splitMethod="bootcv",
                    B=10,
                  data=Melanoma,
                  maxtime=2500)


cv.cindex <- cindex(list(AbsRisk=arr.fit,CauseSpecCox=cox.fit,FineGray=fg.fit,LogisticRisk=lrr.fit),
                  formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,
                  cens.model="cox",
                   splitMethod="bootcv",
                    B=10,
                  data=Melanoma,
                  maxtime=2500)


## graphical model checking
tarr.fit <- ARR(Hist(time,status)~sex+epicel+strata(ulcer)+age+strata(logthick),data=Melanoma,cause=1,
                    cens.model="cox",
                    cens.formula=~sex+epicel+ulcer+age+logthick)
plotEffects(tarr.fit,formula=~ulcer,xlim=c(1000,3000),ylim=c(-3,3))
plotEffects(tarr.fit,formula=~logthick,xlim=c(1000,3000),ylim=c(-3,3))
