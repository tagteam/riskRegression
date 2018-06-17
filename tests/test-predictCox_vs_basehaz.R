### test-predictCox-basehaz.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (10:38) 
## Version: 
## last-updated: jun 16 2018 (11:30) 
##           By: Brice Ozenne
##     Update #: 49
#----------------------------------------------------------------------
## 
### Commentary: 
##
## Check the estimation of the baseline hazard:
##   - compare the output of riskRegression:::predictCox with survival:::basehaz.
##     In particular check that they output the same results when centered = TRUE and centered = FALSE.
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Settings
library(riskRegression)
library(testthat)
library(rms)
library(survival)

context("[predictCox] baseline hazard")
## * [predictCox] baseline hazard (no strata)
cat("[predictCox] Estimation of the baseline hazard (no strata) \n")

## ** Data
data(Melanoma)

## ** Model
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)

## ** Compare to survival::basehaz
test_that("baseline hazard (no strata): compare to survival::basehaz",{
  ## vs basehaz
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
               survival::basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
               survival::basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumhazard, 
               survival::basehaz(fit.cph)$hazard, tolerance = 1e-8)

  ## consistency cph coxph
  ## possible differences due to different fit - coef(fit.coxph)-coef(fit.cph)
  expect_equal(predictCox(fit.cph),
               predictCox(fit.coxph, centered = TRUE), 
               tolerance = 100*max(abs(coef(fit.coxph)-coef(fit.cph))))
})

## ** Number of events
test_that("baseline hazard (no strata): number of events",{

    ## find all unique event times
    GS.alltime <- sort(unique(Melanoma$time))
    
    RR.coxph <- predictCox(fit.coxph)
    expect_equal(RR.coxph$times, GS.alltime)

    RR.cph <- predictCox(fit.cph)
    expect_equal(RR.cph$times, GS.alltime)
})

## * [predictCox] baseline hazard (strata)
cat("[predictCox] Estimation of the baseline hazard (strata) \n")

## ** Data
data(Melanoma)

## ** Model
fitS.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
fitS.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)

## ** Compare to survival::basehaz
test_that("baseline hazard (strata): compare to survival::basehaz",{

  ## vs basehaz
  expect_equal(predictCox(fitS.coxph, centered = FALSE)$cumhazard, 
               basehaz(fitS.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fitS.coxph, centered = TRUE)$cumhazard, 
               basehaz(fitS.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fitS.cph)$cumhazard, 
               basehaz(fitS.cph)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fitS.coxph, keep.strata = TRUE)$strata, 
               basehaz(fitS.coxph)$strata)
  expect_equal(predictCox(fitS.cph, keep.strata = TRUE)$strata, 
               basehaz(fitS.cph)$strata)

  ## consistency cph coxph
  ## different ordering of the strata
  e.coxph <- as.data.table(predictCox(fitS.coxph))
  e.cph <- as.data.table(predictCox(fitS.cph))
  levels(e.coxph$strata)
  levels(e.cph$strata)
})

## ** Number of events
test_that("baseline hazard (strata): number of events",{

    GS.alltime <- tapply(Melanoma$time, interaction(Melanoma$ici,Melanoma$invasion),
                         function(x){sort(unique(x))})

    RR.coxph <- predictCox(fitS.coxph)
    test.alltime <- tapply(RR.coxph$times, RR.coxph$strata, function(x){x})

    expect_equal(unname(GS.alltime),unname(test.alltime))


    GS.alltime <- tapply(Melanoma$time, interaction(Melanoma$invasion,Melanoma$ici),
                         function(x){sort(unique(x))})

    RR.cph <- predictCox(fitS.cph)
    test.alltime <- tapply(RR.cph$times, RR.cph$strata, function(x){x})

    expect_equal(unname(GS.alltime),unname(test.alltime))
})

## * [predictCox] baseline hazard with time varying covariates (no strata)
cat("[predictCox] Estimation of the baseline hazard (time varying cov, no strata) \n")
## ** Data
## example from help(coxph)
dt.TV <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 

## ** Model
fit.coxphTV <- coxph(Surv(start, stop, event) ~ x, data = dt.TV, x = TRUE, y = TRUE)
fit.cphTV <- cph(Surv(start, stop, event) ~ x, data = dt.TV, x = TRUE, y = TRUE)


## ** Compare to survival::basehaz
test_that("baseline hazard (no strata, time varying): compare to survival::basehaz",{

    expect_equal(suppressWarnings(predictCox(fit.coxphTV, centered = FALSE)$cumhazard), 
                 basehaz(fit.coxphTV, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fit.coxphTV, centered = TRUE)$cumhazard), 
                 basehaz(fit.coxphTV, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fit.cphTV)$cumhazard), 
                 basehaz(fit.cphTV)$hazard, tolerance = 1e-8)
  
})

## * [predictCox] baseline hazard with time varying covariates (strata)
cat("[predictCox] Estimation of the baseline hazard (time varying cov, strata) \n")
## ** Data
set.seed(10)
dtS.TV <- rbind(cbind(as.data.table(dt.TV),S = 1),
                      cbind(as.data.table(dt.TV),S = 2))
dtS.TV[, randomS := rbinom(.N,size = 1, prob = 1/2)]

## ** Model
fitS1.coxphTV <- coxph(Surv(start, stop, event) ~ strata(S) + x, data = dtS.TV, x = TRUE, y = TRUE)
fitS1.cphTV <- cph(Surv(start, stop, event) ~ strat(S) + x, data = dtS.TV, x = TRUE, y = TRUE)

fitS2.coxphTV <- coxph(Surv(start, stop, event) ~ strata(randomS) + x, data = dtS.TV, x = TRUE, y = TRUE)
fitS2.cphTV <- cph(Surv(start, stop, event) ~ strat(randomS) + x, data = dtS.TV, x = TRUE, y = TRUE)


## **  Compare to survival::basehaz
test_that("baseline hazard (strata, time varying): compare to survival::basehaz",{

    ## strata defined by S
    expect_equal(suppressWarnings(predictCox(fitS1.coxphTV, centered = FALSE)$cumhazard), 
                 basehaz(fitS1.coxphTV, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fitS1.coxphTV, centered = TRUE)$cumhazard), 
                 basehaz(fitS1.coxphTV, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fitS1.cphTV)$cumhazard), 
                 basehaz(fitS1.cphTV)$hazard, tolerance = 1e-8)

    ## strata defined by randomS
    expect_equal(suppressWarnings(predictCox(fitS2.coxphTV, centered = FALSE)$cumhazard), 
                 basehaz(fitS2.coxphTV, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fitS2.coxphTV, centered = TRUE)$cumhazard), 
                 basehaz(fitS2.coxphTV, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(suppressWarnings(predictCox(fitS2.cphTV)$cumhazard), 
                 basehaz(fitS2.cphTV)$hazard, tolerance = 1e-8)

})



#----------------------------------------------------------------------
### test-predictCox-basehaz.R ends here
