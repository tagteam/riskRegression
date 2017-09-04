### test-predictCox-basehaz.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (10:38) 
## Version: 
## last-updated: sep  4 2017 (10:44) 
##           By: Brice Ozenne
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
##
## Compare the output of riskRegression:::predictCox with survival:::basehaz.
## In particular check that they output the same results when centered = TRUE and centered = FALSE.
## Also perform additional tests on the other argument of predictCox.
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(testthat)
library(rms)
library(survival)

# {{{ 1- [predictCox] Check baseline hazard
cat("Estimation of the baseline hazard \n")
data(Melanoma)

test_that("baseline hazard - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
               basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumhazard, 
               basehaz(fit.cph)$hazard, tolerance = 1e-8)
  
  ## possible differences due to different fit - coef(fit.coxph)-coef(fit.cph)
  expect_equal(predictCox(fit.cph),
               predictCox(fit.coxph, centered = TRUE), 
               tolerance = 100*max(abs(coef(fit.coxph)-coef(fit.cph))))
})

#### strata
test_that("baseline hazard (strata) - order of the results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, keep.strata = TRUE)$strata, 
               basehaz(fit.coxph)$strata)
  expect_equal(predictCox(fit.cph, keep.strata = TRUE)$strata, 
               basehaz(fit.cph)$strata)
  # expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), 
  #              as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata))
})

test_that("baseline hazard (strata) - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumhazard, 
               basehaz(fit.cph)$hazard)
  
})

# }}}

# {{{ 2- [predictCox] Check format of the output
cat("Format of the output \n")
data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

## no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard - correct number of events",{
    # c("time","hazard","cumhazard","survival") remove lastEventTime from pfit
    # time hazard cumhazard survival should have length equals to the number of eventtimes (including censored events)
    # this is not true for lastEventTime which is has length the number of strata
    pfit.coxph <- predictCox(fit.coxph, type = c("hazard","cumhazard","survival"), keep.times = TRUE)[c("time","hazard","cumhazard","survival")]
    lengthRes <- unlist(lapply(pfit.coxph, length))
    expect_equal(unname(lengthRes), rep(length(unique(fit.coxph$y[,"time"])), 4))
    pfit.cph <- predictCox(fit.cph, type = c("hazard","cumhazard","survival"), keep.times = TRUE)[c("time","hazard","cumhazard","survival")]
    lengthRes <- unlist(lapply(pfit.cph, length))
    expect_equal(unname(lengthRes), rep(length(unique(fit.cph$y[,"time"])), 4))
})

## strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard (strata) - order of the results",{
  expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata),
               as.numeric(basehaz(fit.coxph)$strata))
  expect_equal(as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata),
               as.numeric(basehaz(fit.cph)$strata))
})

test_that("baseline hazard (strata) - correct number of events",{
    # c("time","hazard","cumhazard","survival", "strata") remove lastEventTime from pfit
    # time hazard cumhazard survival and strata should have length equals to the number of eventtimes (including censored events)
    # this is not true for lastEventTime which is has length the number of strata

  strata <- interaction(Melanoma$invasion, Melanoma$ici)
  timePerStrata <- tapply(fit.coxph$y[,"time"],strata, function(x){length(unique(x))})

  pfit.coxph <- predictCox(fit.coxph, type = c("hazard","cumhazard","survival"), keep.times = TRUE, keep.strata = TRUE)[c("time","hazard","cumhazard","survival","strata")]
  lengthRes <- unlist(lapply(pfit.coxph, length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
  pfit.cph <- predictCox(fit.cph, type = c("hazard","cumhazard","survival"), keep.times = TRUE, keep.strata = TRUE)[c("time","hazard","cumhazard","survival","strata")]
  lengthRes <- unlist(lapply(pfit.cph, length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
})



test_that("Prediction with Cox model (strata) - export of strata and times",{
    fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
    fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
    predictTempo <- predictCox(fit.coxph)
    expect_equal(length(predictTempo$strata)>0, TRUE) 
    expect_equal(length(predictTempo$time)>0, TRUE)
    predictTempo <- predictCox(fit.coxph, keep.strata = FALSE) # as.data.table(predictCox(fit.coxph, keep.strata = TRUE))
    expect_equal(length(predictTempo$strata)>0, FALSE) 
    expect_equal(length(predictTempo$time)>0, TRUE)
    predictTempo <- predictCox(fit.coxph, keep.strata = FALSE, keep.times = FALSE)
    expect_equal(length(predictTempo$strata)>0, FALSE)
    expect_equal(length(predictTempo$time)>0, FALSE)
  
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1)
  expect_equal(length(predictTempo$strata)>0, TRUE) 
  expect_equal(length(predictTempo$time)>0, TRUE)
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1, keep.strata = FALSE)
  expect_equal(length(predictTempo$strata)>0, FALSE)
  expect_equal(length(predictTempo$time)>0, TRUE)
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1, keep.strata = FALSE, keep.times = FALSE)
  expect_equal(length(predictTempo$strata)>0, FALSE)
  expect_equal(length(predictTempo$time)>0, FALSE)
})

test_that("Prediction with Cox model (strata) - consistency of hazard/cumhazard/survival",{
  predictTempo <- predictCox(fit.coxph, type = c("hazard","cumhazard","survival"), times = times1, newdata = dataset1)
  expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumhazard,1,diff)), tolerance = 1e-8)
  expect_equal(predictTempo$survival, exp(-predictTempo$cumhazard), tolerance = 1e-8)
})

predictTempo <- predictCox(fit.coxph, type = c("hazard","cumhazard","survival"), times = c(0,times1[1:10]), newdata = dataset1[1:2,])
expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumhazard,1,diff)), tolerance = 1e-8)
expect_equal(predictTempo$survival, exp(-predictTempo$cumhazard), tolerance = 1e-8)


test_that("Prediction with Cox model (strata) - incorrect strata",{
    fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
    dataset1$invasion <- "5616"
    expect_error(predictCox(fit.coxph, times = times1, newdata = dataset1))
})
# }}}


#----------------------------------------------------------------------
### test-predictCox-basehaz.R ends here
