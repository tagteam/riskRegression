### test-phreg.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 28 2017 (09:52) 
## Version: 
## last-updated: feb 28 2017 (15:03) 
##           By: Brice Ozenne
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

test <- FALSE

if(test){
library(testthat)
context("Compatibility with the phreg function from the mets package")
library(riskRegression)
library(pec)
library(rms)
library(survival)
library(prodlim)

# {{{ survival

set.seed(10)
d <- sampleData(5e1, outcome = "survival")[,.(eventtime,event,X1,X2,X3,X4,X6)]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]
d[, entry:= 0]
d[, id := 1:.N]
d[, group := factor(rbinom(.N,2,0.5))]

## no strata
m.phreg <- phreg(Surv(entry,eventtime,event)~X1+X2,data=d)
m.coxph <- coxph(Surv(eventtime,event)~X1+X2,data=d, x = TRUE, y = TRUE)
expect_equal(predictCox(m.phreg),predictCox(m.coxph), tol = 1e-3)

predictCox(m.phreg, iid = TRUE, newdata = d[1:5,], times = 1:6) 
outSE$cumhazard.iid

pred.phreg <- predictCox(m.phreg, newdata = d, time = 1:5, se = TRUE)
pred.coxph <- predictCox(m.coxph, newdata = d, time = 1:5, se = TRUE)
expect_equal(pred.phreg$survival, pred.coxph$survival, tol = 1e-8)

## one strata variable
mS.phreg <- phreg(Surv(entry,eventtime,event)~X1+X2+strata(group)+cluster(id),data=d)
mS.coxph <- coxph(Surv(eventtime,event)~X1+X2+strata(group),data=d, x = TRUE, y = TRUE)
expect_equal(predictCox(mS.phreg),predictCox(mS.coxph), tol = 1e-3)

## several strata variables
mS.phreg <- phreg(Surv(entry,eventtime,event)~X1+X2+strata(group,X3)+cluster(id),data=d)
mS.coxph <- coxph(Surv(eventtime,event)~X1+X2+strata(group)+strata(X3),data=d, x = TRUE, y = TRUE)
expect_equal(predictCox(mS.phreg)[c("hazard","cumhazard","survival")],
             predictCox(mS.coxph)[c("hazard","cumhazard","survival")],
             tol = 1e-3)
# }}}

# {{{ Competing risks
set.seed(10)
d <- sampleData(5e1, outcome = "competing.risks")[,.(time,event,X1,X2,X3,X4,X6)]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]
d[, id := 1:.N]
d[, group := factor(rbinom(.N,2,0.5))]

## no strata
m.phreg <- CSC(Hist(time, event)~X1+X2,data=d,
               fitter = "phreg")
m.coxph <- CSC(Hist(time,event)~X1+X2,data=d,
               fitter = "coxph")
expect_equal(class(m.phreg$models[[1]]),"phreg")
pred.phreg <- predict(m.phreg, newdata = d, times = 1:5, cause = 1, se = TRUE)
pred.coxph <- predict(m.coxph, newdata = d, times = 1:5, cause = 1, se = TRUE)
expect_equal(pred.phreg,pred.coxph, tol = 1e-3)

## one strata variable
mS.phreg <- CSC(Hist(time, event)~X1+X2+strata(group)+cluster(id),data=d,
               fitter = "phreg")
mS.coxph <- CSC(Hist(time,event)~X1+X2+strata(group),data=d,
               fitter = "coxph")
expect_equal(class(mS.phreg$models[[1]]),"phreg")
predS.phreg <- predict(mS.phreg, newdata = d, times = 1:5, cause = 1, se = TRUE)
predS.coxph <- predict(mS.coxph, newdata = d, times = 1:5, cause = 1, se = TRUE)
expect_equal(predS.phreg,predS.coxph, tol = 1e-3)

## several strata variables
mS.phreg <- CSC(Hist(time, event)~X1+X2+strata(group,X3)+cluster(id),data=d,
               fitter = "phreg")
mS.coxph <- CSC(Hist(time,event)~X1+X2+strata(group)+strata(X3),data=d,
               fitter = "coxph")
expect_equal(class(mS.phreg$models[[1]]),"phreg")
predS.phreg <- predict(mS.phreg, newdata = d, times = 1:5, cause = 1, se = TRUE)
predS.coxph <- predict(mS.coxph, newdata = d, times = 1:5, cause = 1, se = TRUE)
expect_equal(predS.phreg,predS.coxph, tol = 1e-3)
# }}}
}

#----------------------------------------------------------------------
### test-phreg.R ends here
