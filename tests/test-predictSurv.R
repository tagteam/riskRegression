### test-predictSurv.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 28 2018 (10:06) 
## Version: 
## Last-Updated: sep 25 2019 (12:19) 
##           By: Brice Ozenne
##     Update #: 29
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

## * function predictSurv
## ** Simulate data
set.seed(10)
d <- sampleData(2e2, outcome = "competing.risk")
d.pred <- d[5:15]
seqTime <- c(0,d[["time"]][5:15],2.45,1e5)

## ** hazard CSC
e.CSC <- CSC(Hist(time, event)~ strata(X1) + strata(X2), data = d, surv.type = "hazard")
e.cox <- coxph(Surv(time, event>0)~ strata(X1) + strata(X2), data = d, x = TRUE, y = TRUE)
## e.cox <- coxph(Surv(time, event>0)~ 1, data = d, x = TRUE, y = TRUE)

test_that("predictSurv (hazard)", {
    ## diag = FALSE 
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime,
                    product.limit = FALSE, iid = TRUE, se = TRUE, average.iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime,
                     iid = TRUE, se = TRUE, average.iid = TRUE)

    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

    testPL <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime,
                    product.limit = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE)
    GSPL <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime,
                     iid = TRUE, se = TRUE, average.iid = TRUE)

    ratioSurv <- test$survival/testPL$survival
    testPL$survival.iid2 <- testPL$survival.iid
    for(iObs in 1:NROW(d)){
        testPL$survival.iid2[,,iObs] <- testPL$survival.iid[,,iObs]*ratioSurv
    }
    expect_equal(GSPL$survival,testPL$survival)
    expect_equal(GSPL$survival.se, testPL$survival.se * ratioSurv)
    expect_equal(GSPL$survival.iid,testPL$survival.iid2)
    expect_equal(t(apply(testPL$survival.iid,2:3,mean)),
                 testPL$survival.average.iid)
})

test_that("predictSurv (hazard,diag)", {
    ## diag = FALSE 
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = d.pred$time,
                    diag = TRUE, product.limit = FALSE, iid = TRUE, se = TRUE, average.iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred, times = d.pred$time,
                     diag = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE)

    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

    testPL <- predict(e.CSC, type = "survival", newdata = d.pred, times = d.pred$time,
                    diag = TRUE, product.limit = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE)
    GSPL <- predictCoxPL(e.cox, newdata = d.pred, times = d.pred$time,
                     diag = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE)

    ratioSurv <- test$survival/testPL$survival
    testPL$survival.iid2 <- testPL$survival.iid
    for(iObs in 1:NROW(d)){
        testPL$survival.iid2[,,iObs] <- testPL$survival.iid[,,iObs]*ratioSurv
    }
    expect_equal(GSPL$survival,testPL$survival)
    expect_equal(GSPL$survival.se, testPL$survival.se * ratioSurv)
    expect_equal(GSPL$survival.iid,testPL$survival.iid2)
    expect_equal(t(apply(testPL$survival.iid,2:3,mean)),
                 testPL$survival.average.iid)
})


## ** survival CSC
expect_equal("[predictCSC] vs. predictCox (no strata) - surv.type=\"survival\"",{
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCox(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)
})

expect_equal("[predictCSC] vs. predictCox (strata) - surv.type=\"survival\"",{
    e.CSC <- CSC(Hist(time, event)~ X6 + strata(X1), data = d, surv.type = "survival")
    
    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCox(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)
    
    ## different strata for each cause
    e.CSC <- CSC(list(Hist(time, event)~ X6 + strata(X1),
                      Hist(time, event)~ X6 + strata(X2)),
                 data = d, surv.type = "survival")
    
    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCox(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(GS$survival,test$survival)
    expect_equal(GS$survival.se,test$survival.se)
    expect_equal(GS$survival.iid,test$survival.iid)
    expect_equal(GS$survival.average.iid,test$survival.average.iid)

})


######################################################################
### test-predictSurv.R ends here
