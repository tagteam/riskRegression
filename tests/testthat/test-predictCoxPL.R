### test-predictCoxPL.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 21 2018 (15:10) 
## Version: 
## Last-Updated: Sep 17 2022 (07:02) 
##           By: Thomas Alexander Gerds
##     Update #: 26
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(testthat)
library(survival)

context("function predictCoxPL")

## * Check against survfit
cat("[predictCoxPL] Check baseline hazard estimation vs survfit \n")
set.seed(10)
d <- sampleData(25, outcome = "survival")

test_that("[predictCoxPL] check against survfit, no strata",{

    fit <- coxph(Surv(time,event)~ 1,
                  data=d, ties="breslow", x = TRUE, y = TRUE)
    
    GS <- survfit(Surv(time,event)~1, data = d)
    test <- predictCoxPL(fit)

    expect_equal(ignore_attr=TRUE,test$survival, GS$surv)
})

test_that("[predictCoxPL] check against survfit, strata",{

    fitS <- coxph(Surv(time,event)~ strata(X2),
                  data=d, ties="breslow", x = TRUE, y = TRUE)
    
    GS <- survfit(Surv(time,event)~X2, data = d)
    test <- predictCoxPL(fitS)

    expect_equal(ignore_attr=TRUE,test$survival, GS$surv)
})

## * Argument diag + iid
cat("[predictCoxPL] Check argument \'diag\' + iid vs predictCox \n")
set.seed(10)
dt <- sampleData(50, outcome = "survival")[,.(time,event,X1,X2,X6)]

test_that("[predictCoxPL] diag no strata", {
    e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE,)

    GS0 <- predictCox(e.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE)
    suppressWarnings(GS <- predictCoxPL(e.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE))
    test <- predictCoxPL(e.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE, diag = TRUE)

    ## check hazard/survival
    expect_equal(ignore_attr=TRUE,dt$time, as.double(test$time))
    expect_equal(ignore_attr=TRUE,diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(ignore_attr=TRUE,diag(GS$survival), as.double(test$survival))

    ## check iid for survival
    expect_equal(ignore_attr=TRUE,GS0$survival.iid, GS$survival.iid) ## same predictCox / predictCoxPL
    expect_equal(ignore_attr=TRUE,GS$survival.iid[,1,1], test$survival.iid[,1,1]) ## same when using diag
    expect_equal(ignore_attr=TRUE,t(apply(GS$survival.iid,1,diag)), test$survival.iid[,1,])
    
    ## check average.iid for survival
    expect_equal(ignore_attr=TRUE,GS0$survival.average.iid, GS$survival.average.iid) ## same predictCox / predictCoxPL
    expect_equal(ignore_attr=TRUE,rowMeans(t(apply(GS$survival.iid,1,diag))), test$survival.average.iid[,1])
})

test_that("[predictCoxPL] diag strata", {
    eS.coxph <- coxph(Surv(time, event) ~ strata(X1) + X6, data = dt, y = TRUE, x = TRUE)

    GS0 <- predictCox(eS.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE)
    suppressWarnings(GS <- predictCoxPL(eS.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE))
    test <- predictCoxPL(eS.coxph, newdata = dt, times = dt$time, iid = TRUE, average.iid = TRUE, se = FALSE, diag = TRUE)

    ## check hazard/survival
    expect_equal(ignore_attr=TRUE,dt$time, as.double(test$time))
    expect_equal(ignore_attr=TRUE,diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(ignore_attr=TRUE,diag(GS$survival), as.double(test$survival))

    ## check iid for survival
    expect_equal(ignore_attr=TRUE,GS0$survival.iid, GS$survival.iid) ## same predictCox / predictCoxPL
    expect_equal(ignore_attr=TRUE,GS$survival.iid[,1,1], test$survival.iid[,1,1]) ## same when using diag
    expect_equal(ignore_attr=TRUE,t(apply(GS$survival.iid,1,diag)), test$survival.iid[,1,])

    ## check average.iid for survival
    expect_equal(ignore_attr=TRUE,GS0$survival.average.iid, GS$survival.average.iid) ## same predictCox / predictCoxPL
    expect_equal(ignore_attr=TRUE,rowMeans(t(apply(GS$survival.iid,1,diag))), test$survival.average.iid[,1])
})

######################################################################
### test-predictCoxPL.R ends here
