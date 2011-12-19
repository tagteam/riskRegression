# {{{ load libraries and data
library(riskRegression)
## install from Rforge
## install.packages("riskRegression", repos = "http://R-Forge.R-project.org")
library(pec)    ## to estimate the prediction performance
## install from Rforge
## install.packages("pec", repos = "http://R-Forge.R-project.org")
## may need most recent version of prodlim 
## install.packages("prodlim", repos = "http://R-Forge.R-project.org")
library(cmprsk) ## for comparison with the Fine-Gray model
data(Melanoma)
help(Melanoma)
# }}}

# {{{ fitting different models
## Absolute risk regression
arr.fit <- ARR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(arr.fit)

## Logistic risk regression
lrr.fit <- LRR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(lrr.fit)

## Fine-Gray regression model
fg.fit <- FGR(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma,cause=1)
print(fg.fit)

## combined cause-specific Cox models
cox.fit <- CSC(Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,data=Melanoma)
print(cox.fit)
# }}}

# {{{ 

# }}}

# {{{ model checking for the ARR

## model with time-varying effect for invasion, ulcer and age
times <- sort(unique(Melanoma$time[Melanoma$time>800]))
tarr.fit <- ARR(Hist(time,status)~timevar(invasion)+sex+epicel+timevar(ulcer)+timevar(age)+thick,data=Melanoma,cause=1,times=times)
## comp.fit <- comp.risk(Surv(time,status!=0)~invasion+const(sex)+const(epicel)+ulcer+age+const(thick),data=Melanoma,cause=Melanoma$status,causeS=1,model="rcif",time.pow=0)
## comp.fit$var.cum
## comp.fit$cum
## plot(comp.fit,specific.comp=3)
## plot(comp.fit,score=1,xlim=c(800,3000))

plotCoef(tarr.fit,~invasion,xlim=c(800,3000),ylim=c(-4,4))

plotEffects(tarr.fit,~invasion,xlim=c(800,3000),ylim=c(-4,4))

plotEffects(tarr.fit,~ulcer,xlim=c(800,3000))
plotEffects(tarr.fit,~age,ylim=c(-0.1,0.1),xlim=c(800,3000))

tarr.fit2 <- ARR(Hist(time,status)~invasion+sex+epicel+timevar(ulcer)+age+thick,data=Melanoma,cause=1,times=times)
plotEffects(tarr.fit2,~ulcer,xlim=c(800,3000))
# }}}

# {{{ Show individual predictions

## combined cause-specific Cox models
plotPredictEventProb(cox.fit,
                     newdata=Melanoma,
                     col=as.numeric(Melanoma$ulcer))
## absolute risk regression model
plotPredictEventProb(arr.fit,
                     newdata=Melanoma,
                     col=as.numeric(Melanoma$ulcer))
## logistic risk regression model
plotPredictEventProb(lrr.fit,
                     newdata=Melanoma,
                     col=as.numeric(Melanoma$ulcer))
## Fine-Gray model
plotPredictEventProb(fg.fit,
                     newdata=Melanoma,
                     col=as.numeric(Melanoma$ulcer))

# }}}

# {{{ Comparison of predictions
# e.g. ARR versus combined Cox regression
# after 5 years = 1826.25 days
plot(predictEventProb(arr.fit,times=1826.25,newdata=Melanoma),
     predictEventProb(cox.fit,times=1826.25,newdata=Melanoma),
     main="Predicted probabilities of death due to malignant melanoma\nafter five years",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Absolute risk model",
     ylab="Combined cause-specific Cox models",
     axes=FALSE)
prodlim:::PercentAxis(1,at=seq(0,1,.25))
prodlim:::PercentAxis(2,at=seq(0,1,.25))
abline(a=0,b=1)
# }}}

# {{{ prediction error

## compare alternative link functions
perr.link <- pec(list(AbsRisk=arr.fit,
                       CSC=cox.fit,
                       FineGray=fg.fit,
                       LogisticRisk=lrr.fit),
                  formula=Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,
                  data=Melanoma,
                  maxtime=3000)
plot(perr.link)

## compare absrisk model with time-constant effects
## to absrisk model where ulcer has a time-varying effect
perr.tv <- pec(list("AbsRisk"=arr.fit,
                 "AbsRisk: effect Ulcer timevar"=tarr.fit2),
            formula=Hist(time,status)~invasion+sex+epicel+ulcer+age+thick,
            data=Melanoma,
            maxtime=3000)
plot(perr.tv)

# }}}
