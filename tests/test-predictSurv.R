### test-predictSurv.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 28 2018 (10:06) 
## Version: 
## Last-Updated: sep  8 2019 (23:35) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(rms)
library(survival)
library(testthat)
context("predictSurv checks")

## * Simulate data
set.seed(10)
d <- sampleData(2e2, outcome = "competing.risk")
d.pred <- d[5:15]
seqTime <- c(0,d[["time"]][5:15],2.45,1e5)

## * hazard case
e.CSC <- CSC(Hist(time, event)~ strata(X1) + strata(X2), data = d, surv.type = "hazard")
e.cox <- coxph(Surv(time, event>0)~ strata(X1) + strata(X2), data = d, x = TRUE, y = TRUE)
## e.cox <- coxph(Surv(time, event>0)~ 1, data = d, x = TRUE, y = TRUE)

test_that("predictSurv (hazard)", {
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = FALSE, iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime, iid = TRUE)
    expect_equal(GS$survival, test$survival)
    expect_equal(GS$survival.iid, test$survival.iid)

    GS$survival.iid[,2,]
    test$survival.iid[,2,]

    n.tt <- 2
    test <- predict(e.CSC, type = "survival", newdata = d.pred[1:n.tt], times = seqTime[1:n.tt], product.limit = FALSE, iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred[1:n.tt], times = seqTime[1:n.tt], iid = TRUE)

    test$survival.iid/GS$cumhazard.iid
    GS$survival.iid/GS$cumhazard.iid

    expect_equal(GS$survival, test$survival)
    expect_equal(GS$survival.iid, test$survival.iid)
    expect_equal(GS$cumhazard.iid, test$survival.iid)

    table(GS$cumhazard.iid[,1,1:10]/GS$survival.iid[,1,1:10])
    
    tcrossprod(test$survival.iid[,1,])
    tcrossprod(GS$survival.iid[,1,])

    tcrossprod(test$survival.iid[,2,])
    tcrossprod(GS$survival.iid[,2,])


    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = TRUE, iid = TRUE)
    GS <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime)$survival
    expect_equal(GS, test)
})

## * survival case
e.CSC <- CSC(Hist(time, event)~ X1 + X2, data = d, surv.type = "survival")
e.cox <- coxph(Surv(time, event>0)~ X1 + X2, data = d, x = TRUE, y = TRUE)

test_that("predictSurv (survival)", {
    test <- predictSurv(e.CSC, newdata = d.pred, times = seqTime, product.limit = FALSE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime)$survival
    expect_equal(GS, test)

    test <- predictSurv(e.CSC, newdata = d.pred, times = seqTime, product.limit = TRUE)
    GS <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime)$survival
    expect_equal(GS, test)
})

######################################################################
### test-predictSurv.R ends here
