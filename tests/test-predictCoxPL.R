### test-predictCoxPL.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 21 2018 (15:10) 
## Version: 
## Last-Updated: jun 27 2018 (15:53) 
##           By: Brice Ozenne
##     Update #: 7
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

context("[predictCoxPL] Computation of the baseline hazard - comparison to survival::survfit")

d <- sampleData(1.1e2, outcome = "survival")

## * Check against survfit
cat("[predictCoxPL] Check against survfit \n")
test_that("[predictCoxPL] check against survfit, no strata",{

    fit <- coxph(Surv(time,event)~ 1,
                  data=d, ties="breslow", x = TRUE, y = TRUE)
    
    GS <- survfit(Surv(time,event)~1, data = d)
    test <- predictCoxPL(fit)

    expect_equal(test$survival, GS$surv)
})

test_that("[predictCoxPL] check against survfit, strata",{

    fitS <- coxph(Surv(time,event)~ strata(X2),
                  data=d, ties="breslow", x = TRUE, y = TRUE)
    
    GS <- survfit(Surv(time,event)~X2, data = d)
    test <- predictCoxPL(fitS)

    expect_equal(test$survival, GS$surv)
})

## * Diag argument
cat("[predictCoxPL] diag argument \n")
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]

test_that("[predictCoxPL] diag no strata", {
    e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCoxPL(e.coxph, newdata = dt, times = dt$time, se = FALSE)
    test <- predictCoxPL(e.coxph, newdata = dt, times = dt$time, se = FALSE, diag = TRUE)
    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))
})

test_that("[predictCoxPL] diag strata", {
    eS.coxph <- coxph(Surv(time, event) ~ strata(X1) + X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCoxPL(eS.coxph, newdata = dt, times = dt$time, se = FALSE)
    test <- predictCoxPL(eS.coxph, newdata = dt, times = dt$time, se = FALSE, diag = TRUE)

    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))
})

######################################################################
### test-predictCoxPL.R ends here
