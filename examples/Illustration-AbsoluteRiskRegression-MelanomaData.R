# {{{ load libraries and data
library(riskRegression)
## install from Rforge
## install.packages("riskRegression", repos = "http://R-Forge.R-project.org")
library(pec)    ## to estimate the prediction performance
## install from Rforge
## install.packages("pec", repos = "http://R-Forge.R-project.org")
## may need most recent version of prodlim 
## install.packages("prodlim", repos = "http://R-Forge.R-project.org")
library(cmprsk) ## for comparison with Fine-Gray model
data(Melanoma)
# }}}

# {{{ a single categorical factor

## Stratified Aalen-Johansen 
np.inv <- prodlim(Hist(time,status)~invasion,data=Melanoma)

## Absolute risk regression
arr.inv <- ARR(Hist(time,status)~invasion,data=Melanoma,cause=1)
print(arr.inv)

## Logistic risk regression
lrr.inv <- LRR(Hist(time,status)~invasion,data=Melanoma,cause=1)
print(lrr.inv)

## Fine-Gray regression model
fg.inv <- FGR(Hist(time,status)~invasion,data=Melanoma,cause=1)
print(fg.inv)

## combined cause-specific Cox models
cox.inv <- cumincCox(Hist(time,status)~invasion,data=Melanoma)
print(cox.inv)

## comparison of predictions
## e.g. ARR versus combined Cox regression
## after 5 years = 1826.25 days
newd <- data.frame(invasion=levels(Melanoma$invasion))
p.arr <- predictEventProb(arr.inv,times=1826.25,newdata= newd)
p.lrr <- predictEventProb(lrr.inv,times=1826.25,newdata= newd)
p.fg <- predictEventProb(fg.inv,times=1826.25,newdata= newd)
p.cox <- predictEventProb(cox.inv,times=1826.25,newdata= newd)
x <- data.frame(newd,
           ARR=round(100*p.arr,1),
           LRR=round(100*p.lrr,1),
           FGR=round(100*p.fg,1),
           CoxR=round(100*p.cox,1))
cat("\nPredicted probabilities of death\n due to malignant melanoma after 5 years\n\n")
print(x)



## prediction error
perr <- pec(list(Strat.AJ=np.inv,AbsRisk=arr.inv,LogisticRisk=lrr.inv,FineGray=fg.inv),formula=Hist(time,status)~invasion,data=Melanoma,maxtime=3000)
plot(perr)
## conclusion: all models have similar prediction performance

# }}}

# {{{ combining several categorical factors

## Absolute risk regression
arr.strat <- ARR(Hist(time,status)~invasion+sex+epicel+ulcer,data=Melanoma,cause=1)
print(arr.strat)

## Logistic risk regression
lrr.strat <- LRR(Hist(time,status)~invasion+sex+epicel+ulcer,data=Melanoma,cause=1)
print(lrr.strat)

## Fine-Gray regression model
fg.strat <- FGR(Hist(time,status)~invasion+sex+epicel+ulcer,data=Melanoma,cause=1)
print(fg.strat)

## combined cause-specific Cox models
cox.strat <- cumincCox(Hist(time,status)~invasion+sex+epicel+ulcer,data=Melanoma)
print(cox.strat)

## comparison of predictions
## e.g. ARR versus combined Cox regression
## after 5 years = 1826.25 days
plot(predictEventProb(arr.strat,times=1826.25,newdata=Melanoma),
     predictEventProb(cox.strat,times=1826.25,newdata=Melanoma),
     main="Predicted probabilities of death due to malignant melanoma\nafter five years",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Absolute risk model",ylab="Combined cause-specific Cox models",
     axes=FALSE)
prodlim:::PercentAxis(1,at=seq(0,1,.25))
prodlim:::PercentAxis(2,at=seq(0,1,.25))
abline(a=0,b=1)

## prediction error
perr.strat <- pec(list(AbsRisk=arr.strat,LogisticRisk=lrr.strat,FineGray=fg.strat),formula=Hist(time,status)~invasion+sex+epicel+ulcer,data=Melanoma,maxtime=3000)
summary(perr.strat,times=seq(0,3000,365.25))
plot(perr.strat)

# }}}

# {{{ a single continuous variable

## Absolute risk regression
## 
arr.thick <- ARR(Hist(time,status)~thick,data=Melanoma,cause=1)
arr.thick.timevar <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
arr.thick.timepower <- ARR(Hist(time,status)~const(thick,power=1),data=Melanoma,cause=1)
print(arr.thick)
print(arr.thick.timevar)
print(arr.thick.timepower)

## Logistic risk regression
lrr.thick <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
print(lrr.thick)

## Fine-Gray regression model
fg.thick <- FGR(Hist(time,status)~thick,data=Melanoma,cause=1)
print(fg.thick)

## combined cause-specific Cox models
cox.thick <- cumincCox(Hist(time,status)~thick,data=Melanoma)
print(cox.thick)

## prediction error
## different versions of 
perr.thick.arr <- pec(list(ARR.const=arr.thick,ARR.timevar=arr.thick.timevar,ARR.timepower=arr.thick.timepower),formula=Hist(time,status)~thick,data=Melanoma,maxtime=3000)
plot(perr.thick.arr)
## conclusion: model with thick * time does not fit very well

perr.thick <- pec(list(AbsRisk=arr.thick,LogisticRisk=lrr.thick,FineGray=fg.thick),formula=Hist(time,status)~thick,data=Melanoma,maxtime=3000)
plot(perr.thick)
## conclusion: all models have similar prediction performance

# }}}

# {{{ mixture of categorical and continuous variables 

## Absolute risk regression
arr.multi <- ARR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(arr.multi)

## Logistic risk regression
lrr.multi <- LRR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(lrr.multi)

## Fine-Gray regression model
fg.multi <- FGR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(fg.multi)

## combined cause-specific Cox models
cox.multi <- cumincCox(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma)
print(cox.multi)


## comparison of predictions
## e.g. ARR versus combined Cox regression
## after 5 years = 1826.25 days
plot(predictEventProb(arr.multi,times=1826.25,newdata=Melanoma),
     predictEventProb(cox.multi,times=1826.25,newdata=Melanoma),
     main="Predicted probabilities of death due to malignant melanoma\nafter five years",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Absolute risk model",ylab="Combined cause-specific Cox models",
     axes=FALSE)
prodlim:::PercentAxis(1,at=seq(0,1,.25))
prodlim:::PercentAxis(2,at=seq(0,1,.25))
abline(a=0,b=1)


## prediction error

perr.multi <- pec(list(AbsRisk=arr.multi,LogisticRisk=lrr.multi,FineGray=fg.multi),formula=Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,maxtime=3000)
plot(perr.multi)

# }}}
