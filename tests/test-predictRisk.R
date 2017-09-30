### test-predictRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 10 2017 (08:56) 
## Version: 
## Last-Updated: Sep 30 2017 (18:32) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(randomForestSRC)

# {{{ missing data
test_that("Additional arguments: example with imputation of missing data", {
    data(pbc,package="survival")
    set.seed(10)
    forest <- rfsrc(Surv(time,status)~chol+age+sex,data=pbc,ntree=10,nsplit=10)
    ## no specification
    expect_error(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km"))
    ## correct specification
    expect_output(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km",predictRisk.args=list("rfsrc"=list(na.action="na.impute")))))
    ## wrong specification
    expect_error(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km",predictRisk.args=list("randomForestSRC"=list(na.action="na.impute")))))
})
# }}}

# {{{ predictCox/predictCSC - to be reorganized (there is no clear test in this section)
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
method.ties <- "efron"
cause <- 1

n <- 3
set.seed(3)
dn <- SimCompRisk(n)
dn$time <- round(dn$time,2)
dn$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
CSC.h3 <- CSC(Hist(time,event) ~ X1 + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph")
CSC.h1 <- CSC(Hist(time,event) ~ strat(X1) + X3 + X2, data = df.S, ties = method.ties, fitter = "cph")
CSC.h <- CSC(Hist(time,event) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph")
CSC.s <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties, fitter = "coxph")
predictRisk(CSC.h1, newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.h3, newdata = dn, times = c(5,10,15,20), cause = cause)

CSC.h0 <- CSC(Hist(time,event) ~ X1 + X3 + X2, data = df.S, ties = method.ties, fitter = "cph")
predictRisk(CSC.h0, newdata = dn, times = c(5,10,15,20), cause = cause)

predictRisk(CSC.h1, newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s, newdata = dn, times = c(5,10,15,20), cause = cause)


df.S[df.S$time==6.55,c("time","event")]
predictCox(CSC.h$models[[1]],newdata = dn[1,],times=c(2.29,6.55),type="hazard")
predictCox(CSC.h$models[[2]],newdata = dn[1,],times=6.55,type="hazard")

predictCox(CSC.h$models[["Cause 1"]],newdata = dn[1,],times=CSC.h$eventTimes,type="hazard")$hazard

test_that("Prediction with CSC - categorical cause",{
predictRisk(CSC.h, newdata = dn[1,], times = c(2), cause = "1")
})

predictRisk(CSC.h, newdata = dn, times = c(1,2,3.24,3.25,3.26,5,10,15,20), cause = cause)
predictRisk(CSC.s, newdata = dn, times = c(1,5,10,15,20), cause = cause)

predictCox(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20))

predictRisk(CSC.h$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)

predictRisk(CSC.h$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)
# }}}

######################################################################
### test-predictRisk.R ends here
