### test-predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (10:38) 
## Version: 
## last-updated: okt  7 2019 (16:27) 
##           By: Brice Ozenne
##     Update #: 129
#----------------------------------------------------------------------
## 
### Commentary: 
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
library(data.table)
library(timereg); nsim.band <- 500;
context("function predictCox")


## * [predictCox] Baseline hazard (no strata)
cat("[predictCox] Estimation of the baseline hazard (no strata) \n")

## ** Data
data(Melanoma, package = "riskRegression")

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

## * [predictCox] Baseline hazard (strata)
cat("[predictCox] Estimation of the baseline hazard (strata) \n")

## ** Data
data(Melanoma, package = "riskRegression")

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

## * [predictCox] Baseline hazard with time varying covariates (no strata)
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

## * [predictCox] Baseline hazard with time varying covariates (strata)
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


## ** Compare to survival::basehaz
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

## * [predictCox] Predictions,se,band,average.iid (no strata, continuous variables)
cat("[predictCox] Predictions (no strata, continuous) \n")

## ** Data
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]
dt[,X1:=as.numeric(as.character(X1))]
dt[,X2:=as.numeric(as.character(X2))]
dt[ , X16 := X1*X6]

## sorted dataset
dt.sort <- copy(dt)
setkeyv(dt.sort,c("time")) 

## ** Model
e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)
e.cph <- cph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)
e.coxph_sort <- coxph(Surv(time, event) ~ X1*X6, data = dt.sort, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X6) + prop(X1*X6),
                       data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

## ** Consistency between hazard/cumhazard/survival

test_that("[predictCox] - consistency of hazard/cumhazard/survival",{
  predRR <- predictCox(e.coxph, type = c("hazard","cumhazard","survival"), times = sort(dt$time), newdata = dt)
  expect_equal(predRR$hazard[,-1], t(apply(predRR$cumhazard,1,diff)), tolerance = 1e-8)
  expect_equal(predRR$survival, exp(-predRR$cumhazard), tolerance = 1e-8)
})

## ** One time
## *** Extract information
predGS <- predict(e.timereg, newdata = dt, times = 10)
predRR1 <- predictCox(e.coxph, newdata = dt, times = 10,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
predRR2 <- predictCox(e.coxph_sort, newdata = dt, times = 10,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)

predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10,
                        nsim.band = nsim.band)
predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10,
                          nsim.band = nsim.band)

## *** Test vs. cph
test_that("[predictCox] compare survival and survival.se coxph/cph (1 fixed time)",{
    res.cph <- predictCox(e.cph, newdata = dt, times = 10,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    expect_equal(as.data.table(predRR1),
                 as.data.table(res.cph),
                 tol = 1e-3)

})

## *** Test vs. timereg
test_that("[predictCox] compare survival and survival.se to timereg (1 fixed time)",{
    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))
})

## *** Test vs. known result
test_that("[confint.predictCox] compare to known values (1 fixed time, no transformation)",{

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.19073, 0.016148), 
                     "survival.se" = c(0.088812, 0.023259), 
                     "survival.lower" = c(0.016662, 0), 
                     "survival.upper" = c(0.364797, 0.061734), 
                     "survival.quantileBand" = c(2.021964, 2.03706), 
                     "survival.lowerBand" = c(0.011156, 0), 
                     "survival.upperBand" = c(0.370304, 0.063527))

    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)
})

test_that("[confint.predictCox] compare to known values (1 fixed time, log log transformation)",{
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.19073, 0.01615), 
                     "survival.se" = c(0.088812, 0.023259), 
                     "survival.lower" = c(0.05646, 0.00028), 
                     "survival.upper" = c(0.38475, 0.12474), 
                     "survival.quantileBand" = c(2.02196, 2.03706), 
                     "survival.lowerBand" = c(0.05368, 0.00022), 
                     "survival.upperBand" = c(0.39115, 0.13183))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE],
                 GS, tol = 1e-4, scale = 1)
})

## ** At event times

## *** Extract information

vec.time <- sort(dt$time[1:10])
set.seed(10)
predGS <- predict(e.timereg, newdata = dt, times = vec.time)
predRR1 <- predictCox(e.coxph, newdata = dt, times = vec.time,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
predRR2 <- predictCox(e.coxph_sort, newdata = dt, times = vec.time,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)

predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10, nsim.band = nsim.band)
predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10, nsim.band = nsim.band)


## *** Test vs. cph
test_that("[predictCox] compare survival and survival.se coxph/cph (eventtimes)",{
    res.cph <- predictCox(e.cph, newdata = dt, times = vec.time,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    expect_equal(as.data.table(predRR1),
                 as.data.table(res.cph),
                 tol = 1e-3)

})
## *** Test vs. timereg
test_that("[predictCox] compare survival and survival.se to timereg (eventtimes)",{
    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))
})

## *** Test vs. known result
test_that("[confint.predictCox] compare to known values (eventtimes, no transformation)",{

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.170031, 0.170031), 
                     "survival" = c(0.993101, 0.982909), 
                     "survival.se" = c(0.007437, 0.018876), 
                     "survival.lower" = c(0.978526, 0.945913), 
                     "survival.upper" = c(1, 1), 
                     "survival.quantileBand" = c(2.666239, 2.66367), 
                     "survival.lowerBand" = c(0.973273, 0.93263), 
                     "survival.upperBand" = c(1, 1))

    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)
})

test_that("[confint.predictCox] compare to known values (eventtimes, log log transformation)",{
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.17003, 0.17003), 
                     "survival" = c(0.9931, 0.98291), 
                     "survival.se" = c(0.007437, 0.018876), 
                     "survival.lower" = c(0.94395, 0.85811), 
                     "survival.upper" = c(0.99917, 0.99806), 
                     "survival.quantileBand" = c(2.66624, 2.66367), 
                     "survival.lowerBand" = c(0.88353, 0.71525), 
                     "survival.upperBand" = c(0.99961, 0.99911))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE],
                 GS, tol = 1e-4, scale = 1)
})

## ** After last event
test_that("[predictCox] after the last event",{

    predRR1 <- predictCox(e.coxph, newdata = dt, times = 1e8,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    predRR1 <- confint(predRR1, nsim.band = nsim.band)
    
    expect_true(all(is.na(predRR1$survival)))
    expect_true(all(is.na(predRR1$survival.se)))
    expect_true(all(is.na(predRR1$survival.lower)))
    expect_true(all(is.na(predRR1$survival.upper)))
    expect_true(all(is.na(predRR1$survival.lowerBand)))
    expect_true(all(is.na(predRR1$survival.upperBand)))

    vec.time <- max(dt$time) + c(-1e-1,0,1e-1)
    predRR <- predictCox(e.coxph, type = c("hazard","cumhazard","survival"),
                         times = vec.time, newdata = dt)
    M.true <- matrix(c(FALSE, FALSE, TRUE), nrow = NROW(dt), ncol = 3, byrow = TRUE)
    expect_equal(is.na(predRR$hazard), M.true)
    expect_equal(is.na(predRR$cumhazard), M.true)
    expect_equal(is.na(predRR$survival), M.true)

})    

test_that("Prediction - last event censored",{

    dt.lastC <- copy(dt)
    Utimes <- sort(unique(dt$time))
    n.Utimes <- length(Utimes)
    dt.lastC[time==max(time), event := 0]
    
    fit <- coxph(Surv(time, event == 1) ~ X1 + X2 + X6, data = dt.lastC, y = TRUE, x = TRUE)
    predictRR <- predictCox(fit, newdata = dt.lastC, times = tail(Utimes, 2))

    survLast <- predictRR$survival[,2]
    survLastM1 <- predictRR$survival[,1]
    expect_true(all(survLast-survLastM1==0))
    
})

test_that("Prediction - last event death",{

    dt.lastD <- copy(dt)
    Utimes <- sort(unique(dt$time))
    n.Utimes <- length(Utimes)
    dt.lastD[time==max(time), event := 1]
    
    fit <- coxph(Surv(time, event == 1) ~ X1 + X2 + X6, data = dt.lastD, y = TRUE, x = TRUE)
    predictRR <- predictCox(fit, newdata = dt.lastD[1], times = tail(Utimes, 2))

    survLast <- predictRR$survival[,2]
    survLastM1 <- predictRR$survival[,1]
    expect_true(all(survLast<survLastM1))
    
})


## ** Before first event
test_that("[predictCox] before the first event",{

    predRR1 <- predictCox(e.coxph, newdata = dt, times = 1e-8,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    predRR1 <- confint(predRR1, nsim.band = nsim.band)
    
    expect_true(all(predRR1$survival==1))
    expect_true(all(predRR1$cumhazard==0))
    expect_true(all(predRR1$survival.se==0))
    expect_true(all(predRR1$survival.lower==1))
    expect_true(all(predRR1$survival.upper==1))
    expect_true(all(predRR1$survival.lowerBand==1))
    expect_true(all(predRR1$survival.upperBand==1))

})    

## ** Sorted vs. unsorted times

test_that("[predictCox] - sorted vs. unsorted times (no strata)",{
    vec.time <- dt$time
    index.sort <- order(vec.time)
    vec.time_sort <- dt$time[index.sort]
    
    predRR1 <- predictCox(e.coxph, newdata = dt, times = vec.time, se = TRUE)
    predRR2 <- predictCox(e.coxph, newdata = dt, times = vec.time_sort, se = TRUE)

    expect_equal(predRR1$survival[,index.sort],
                 predRR2$survival)

    expect_equal(predRR1$survival.se[,index.sort],
                 predRR2$survival.se)

    expect_equal(predRR1$time[index.sort],
                 predRR2$time)

})

## ** iid.average
test_that("[predictCox] - fast iid average (no strata)",{

    ## simple average
    predRR.av <- predictCox(e.coxph, times = dt$time[1:5], average.iid = TRUE, newdata = dt,
                            type = c("hazard","cumhazard","survival"))
    predRR.GS <- predictCox(e.coxph, times = dt$time[1:5], iid = TRUE, newdata = dt,
                            type = c("hazard","cumhazard","survival"))

    expect_equal(t(apply(predRR.GS$hazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$hazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$cumhazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$survival.iid, MARGIN = 2:3,mean)),
                 predRR.av$survival.average.iid, tolerance = 1e-8)

    ## weighted average
    fT <- TRUE
    attr(fT, "factor") <- list(matrix(5, nrow = NROW(dt), ncol = 5),
                               matrix(1:NROW(dt), nrow = NROW(dt), ncol = 5)
                               )
    
    predRR.av2 <- predictCox(e.coxph, times = sort(dt$time[1:5]), average.iid = fT, newdata = dt,
                             type = c("hazard","cumhazard","survival"))
    calcGS <- function(iid){
        t(apply(iid, MARGIN = 2:3, function(iCol){
            mean(iCol * attr(fT, "factor")[[2]][,1])
        }))
    }

    expect_equal(predRR.av$hazard.average.iid[,order(dt$time[1:5])]*5,
                 predRR.av2$hazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(predRR.av$cumhazard.average.iid[,order(dt$time[1:5])]*5,
                 predRR.av2$cumhazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(predRR.av$survival.average.iid[,order(dt$time[1:5])]*5,
                 predRR.av2$survival.average.iid[[1]],
                 tol = 1e-8)

    expect_equal(calcGS(predRR.GS$hazard.iid)[,order(dt$time[1:5])],
                 predRR.av2$hazard.average.iid[[2]],
                 tolerance = 1e-8)
    expect_equal(calcGS(predRR.GS$cumhazard.iid)[,order(dt$time[1:5])],
                 predRR.av2$cumhazard.average.iid[[2]],
                 tolerance = 1e-8)
    expect_equal(calcGS(predRR.GS$survival.iid)[,order(dt$time[1:5])],
                 predRR.av2$survival.average.iid[[2]],
                 tolerance = 1e-8)
})

## * [predictCox] Predictions,se,band,average.iid (no strata, categorical variable)
cat("[predictCox] Predictions (no strata, categorical) \n")

## ** Data
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]
dt[,X1:=as.numeric(as.character(X1))]
dt[,X2:=as.numeric(as.character(X2))]
dt[ , Xcat2 := as.factor(paste0(X1,X2))]

## sorted dataset
dt.sort <- copy(dt)
setkeyv(dt.sort,c("time")) 

## ** Model
e.coxph <- coxph(Surv(time, event) ~ Xcat2 + X6 , data = dt, y = TRUE, x = TRUE)
e.cph <- cph(Surv(time, event) ~ Xcat2 + X6 , data = dt, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(Xcat2) + prop(X6),
                       data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)


## ** One time and event times
test_that("[predictCox] compare survival and survival.se to timereg/cph (categorical variable)",{
    ## fixed time
    predGS <- predict(e.timereg, newdata = dt, times = 10)
    predRR1 <- predictCox(e.coxph, newdata = dt, times = 10, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    predCPH <- predictCox(e.cph, newdata = dt, times = 10, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predCPH$survival), tol = 1e-3)
    expect_equal(as.double(predRR1$survival.se), as.double(predCPH$survival.se), tol = 1e-3)

    ## event time
    predGS <- predict(e.timereg, newdata = dt, times = sort(dt$time))
    predRR1 <- predictCox(e.coxph, newdata = dt, times = sort(dt$time), se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    predCPH <- predictCox(e.cph, newdata = dt, times = sort(dt$time), se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predCPH$survival), tol = 1e-3)
    expect_equal(as.double(predRR1$survival.se), as.double(predCPH$survival.se), tol = 1e-3)
})



## * [predictCox] Predictions,se,band,average.iid (strata)
cat("[predictCox] Predictions (strata) \n")
## ** Data
set.seed(10)
dtStrata <- copy(dt)
dtStrata[, strata :=  rbinom(n = .N, size = 2, prob = c(1/3,1/2))] # strata
dtStrata.sort <- copy(dtStrata)
setkeyv(dtStrata.sort, c("strata", "time"))

## ** Model
eS.timereg <- cox.aalen(Surv(time, event) ~ strata(strata)-1 + prop(X1) + prop(X6),
                       data = dtStrata, 
                       resample.iid = TRUE, max.timepoint.sim=NULL)

eS.cph <- cph(Surv(time, event) ~ strat(strata) + X1 + X6,
              data = dtStrata, y = TRUE, x = TRUE)

eS.coxph <- coxph(Surv(time, event) ~ strata(strata) + X1 + X6,
                 data = dtStrata, y = TRUE, x = TRUE)

## ** Reject incorrect strata
test_that("[predictCox] - incorrect strata",{
    dt2 <- copy(dt)
    dt2$strata <- "5616"
    expect_error(suppressWarnings(predictCox(eS.coxph, times = dt2$time, newdata = dt2)))
})

## ** 1 fixed time

## *** Extract information
set.seed(10)
predGS <- predict(eS.timereg, newdata = dtStrata, times = 4)
predRR1 <- predictCox(eS.coxph, newdata = dtStrata, times = 4,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)

predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10,
                        nsim.band = nsim.band)
predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10,
                          nsim.band = nsim.band)

## *** Test vs. cph
test_that("[predictCox] compare survival and survival.se coxph/cph (1 fixed time, strata)",{
    set.seed(10)
    res.cph <- predictCox(eS.cph, newdata = dtStrata, times = 4,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    expect_equal(as.data.table(predRR1),
                 as.data.table(res.cph),
                 tol = 1e-3)

})


## *** Test vs. timereg
test_that("[predictCox] compare survival and survival.se to timereg (1 fixed time, strata)",{
    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})


## *** Test vs. known result
test_that("[confint.predictCox] compare to known values (1 fixed time, no transformation, strata)",{

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(4, 4), 
                     "survival" = c(0.780713, 0.519771), 
                     "survival.se" = c(0.082657, 0.12765), 
                     "survival.lower" = c(0.618708, 0.269583), 
                     "survival.upper" = c(0.942717, 0.76996), 
                     "survival.quantileBand" = c(1.855252, 1.944986), 
                     "survival.lowerBand" = c(0.627364, 0.271494), 
                     "survival.upperBand" = c(0.934062, 0.768048))

    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)
})

test_that("[confint.predictCox] compare to known values (1 fixed time, log log transformation, strata)",{
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(4, 4), 
                     "survival" = c(0.78071, 0.51977), 
                     "survival.se" = c(0.082657, 0.12765), 
                     "survival.lower" = c(0.56416, 0.25526), 
                     "survival.upper" = c(0.89848, 0.73082), 
                     "survival.quantileBand" = c(1.85525, 1.94499), 
                     "survival.lowerBand" = c(0.57849, 0.25722), 
                     "survival.upperBand" = c(0.89408, 0.72953))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE],
                 GS, tol = 1e-4, scale = 1)
})

## ** At event times (in one strata)

## *** Extract information

vec.time <- sort(dtStrata$time)[1:10]
set.seed(10)
predGS <- predict(eS.timereg, newdata = dtStrata, times = vec.time)
predRR1 <- predictCox(eS.coxph, newdata = dtStrata, times = vec.time,
                      se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)

predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10,
                        nsim.band = nsim.band)
predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10,
                          nsim.band = nsim.band)

## *** Test vs. cph
test_that("[predictCox] compare survival and survival.se coxph/cph (eventtime, strata)",{
    set.seed(10)
    res.cph <- predictCox(eS.cph, newdata = dtStrata, times = vec.time,
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    expect_equal(as.data.table(predRR1),
                 as.data.table(res.cph),
                 tol = 1e-3)

})

## *** Test vs. timereg
test_that("[predictCox] compare survival and survival.se to timereg (eventtimes, strata)",{
    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})

## *** Test vs. known result
test_that("[confint.predictCox] compare to known values (eventtimes, no transformation, strata)",{

    ## confidence intervals/band
    GS <- data.table("observation" = c(35, 36), 
                     "times" = c(0.877076, 0.877076), 
                     "survival" = c(0.994832, 0.992565), 
                     "survival.se" = c(0.00369, 0.010196), 
                     "survival.lower" = c(0.987599, 0.97258), 
                     "survival.upper" = c(1, 1), 
                     "survival.quantileBand" = c(2.276138, 2.331861), 
                     "survival.lowerBand" = c(0.986432, 0.968788), 
                     "survival.upperBand" = c(1, 1))

    ## butils::object2script(as.data.table(predRR1.none)[135:136,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[135:136,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)
})

test_that("[confint.predictCox] compare to known values (eventtimes, log log transformation, strata)",{
    GS <- data.table("observation" = c(35, 36), 
                     "times" = c(0.87708, 0.87708), 
                     "survival" = c(0.99483, 0.99256), 
                     "survival.se" = c(0.00369, 0.010196), 
                     "survival.lower" = c(0.97914, 0.89511), 
                     "survival.upper" = c(0.99873, 0.9995), 
                     "survival.quantileBand" = c(2.27614, 2.33186), 
                     "survival.lowerBand" = c(0.97391, 0.83119), 
                     "survival.upperBand" = c(0.99898, 0.9997))
    ## butils::object2script(as.data.table(predRR1.loglog)[135:136,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[135:136,names(GS),with=FALSE],
                 GS, tol = 1e-4, scale = 1)
})


## ** after last event
test_that("[predictCox] after the last event (strata)",{
    lastevent <- dtStrata[, max(time), by = "strata"]
    laststrata <- lastevent[V1==max(V1),strata]
    
    predRR1 <- predictCox(eS.coxph, newdata = dtStrata, times = max(lastevent[["V1"]]),
                          se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
    predRR1 <- confint(predRR1, nsim.band = nsim.band)
    
    expect_true(all(is.na(predRR1$survival[dtStrata$strata!=laststrata,])))
    expect_true(all(is.na(predRR1$cumhazard[dtStrata$strata!=laststrata,])))

    expect_true(all(!is.na(predRR1$survival[dtStrata$strata==laststrata,])))
    expect_true(all(!is.na(predRR1$cumhazard[dtStrata$strata==laststrata,])))

    expect_true(all(is.na(predRR1$survival.se[dtStrata$strata!=laststrata,])))
    expect_true(all(is.na(predRR1$survival.lower[dtStrata$strata!=laststrata,])))
    expect_true(all(is.na(predRR1$survival.upper[dtStrata$strata!=laststrata,])))
    expect_true(all(is.na(predRR1$survival.lowerBand[dtStrata$strata!=laststrata,])))
    expect_true(all(is.na(predRR1$survival.upperBand[dtStrata$strata!=laststrata,])))


})    

## ** before first event
test_that("[predictCox] before the first event (strata)",{
    firstevent <- dtStrata[, min(time), by = "strata"]
    firststrata <- firstevent[V1==min(V1),strata]

    predRR1 <- predictCox(eS.coxph, newdata = dtStrata, times = min(firstevent[["V1"]]))
    
    expect_true(all(predRR1$survival[dtStrata$strata!=firststrata,]==1))
    expect_true(all(predRR1$cumhazard[dtStrata$strata!=firststrata,]==0))

    expect_true(all(predRR1$survival[dtStrata$strata==firststrata,]<1))
    expect_true(all(predRR1$cumhazard[dtStrata$strata==firststrata,]>0))

})
## ** iid.average
test_that("[predictCox] - iid average",{
    ## eS.coxph <- coxph(Surv(time, event) ~ strata(strata), data = dtStrata, x = TRUE)
    ## seqTime <- c(0,sort(dtStrata$time)[1:5])
    seqTime <- c(0,dtStrata$time[1:5])
    
    ## simple average
    predRR.av <- predictCox(eS.coxph, times = seqTime, average.iid = TRUE, newdata = dtStrata,
                            type = c("hazard","cumhazard","survival"))
    predRR.GS <- predictCox(eS.coxph, times = seqTime, iid = TRUE, newdata = dtStrata,
                            type = c("hazard","cumhazard","survival"))

    expect_equal(t(apply(predRR.GS$hazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$hazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$cumhazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$survival.iid, MARGIN = 2:3,mean)),
                 predRR.av$survival.average.iid, tolerance = 1e-8)

    ## weighted average
    fT <- TRUE
    attr(fT, "factor") <- list(matrix(5, nrow = NROW(dtStrata), ncol = 6),
                               matrix(1:NROW(dtStrata), nrow = NROW(dtStrata), ncol = 6)
                               )
    
    predRR.av2 <- predictCox(eS.coxph, times = sort(seqTime), average.iid = fT, newdata = dtStrata,
                             type = c("hazard","cumhazard","survival"))
    calcGS <- function(iid){
        t(apply(iid, MARGIN = 2:3, function(iCol){
            mean(iCol * attr(fT, "factor")[[2]][,1])
        }))
    }

    expect_equal(predRR.av$hazard.average.iid[,order(seqTime)]*5,
                 predRR.av2$hazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(predRR.av$cumhazard.average.iid[,order(seqTime)]*5,
                 predRR.av2$cumhazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(predRR.av$survival.average.iid[,order(seqTime)]*5,
                 predRR.av2$survival.average.iid[[1]],
                 tol = 1e-8)

    expect_equal(calcGS(predRR.GS$hazard.iid)[,order(seqTime)],
                 predRR.av2$hazard.average.iid[[2]],
                 tolerance = 1e-8)
    expect_equal(calcGS(predRR.GS$cumhazard.iid)[,order(seqTime)],
                 predRR.av2$cumhazard.average.iid[[2]],
                 tolerance = 1e-8)
    expect_equal(calcGS(predRR.GS$survival.iid)[,order(seqTime)],
                 predRR.av2$survival.average.iid[[2]],
                 tolerance = 1e-8)
})

## * [predictCox] SE/CI check against manual computation
cat("[predictCox] SE/CI check against manual computation \n")
## from confint.predictCox

## ** Data
set.seed(10)
dt <- sampleData(40,outcome="survival") 
 
## ** Model
fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
             data=dt, ties="breslow", x = TRUE, y = TRUE)

fit.pred <- predictCox(fit, newdata=dt[1:3], times=c(3,8), type = "survival",
                       se = TRUE, iid = TRUE, band = TRUE, confint = FALSE)
confint.pred1 <- confint(fit.pred, survival.transform = "none", nsim.band = nsim.band)
confint.pred2 <- confint(fit.pred, survival.transform = "loglog", nsim.band = nsim.band)

## ** standard errors
test_that("[predictCox] consistency  iid se", {
    expect_equal(fit.pred$survival.se[1,],
                 sqrt(rowSums(fit.pred$survival.iid[1,,]^2))
                 )
})

## ** confidence intervals / bands computed on the original scale
test_that("[confint.predictCox] manual computation on ci", {
    ## orignial scale
    expect_equal(confint.pred1$survival.lower,
                 fit.pred$survival + qnorm(0.025) * fit.pred$survival.se)
    expect_equal(as.double(confint.pred1$survival.upper),
                 pmin(1,fit.pred$survival + qnorm(0.975) * fit.pred$survival.se))

    ## loglog scale
    newse <- fit.pred$survival.se/(-fit.pred$survival*log(fit.pred$survival))
    expect_equal(confint.pred2$survival.lower,
                 exp(-exp(log(-log(fit.pred$survival)) + qnorm(0.975) * newse)))
    expect_equal(confint.pred2$survival.upper,
                 exp(-exp(log(-log(fit.pred$survival)) + qnorm(0.025) * newse)))
})

## * [predictCox] Diag argument
cat("[predictCox] Argument \'diag\' \n")
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]

test_that("[predictCox] diag no strata", {
    e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCox(e.coxph, newdata = dt, times = dt$time, se = TRUE, iid = TRUE, average.iid = TRUE)
    test <- predictCox(e.coxph, newdata = dt, times = dt$time,
                       se = TRUE, iid = TRUE, average.iid = TRUE, diag = TRUE)
    test2 <- predictCox(e.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE)

    ## estimates
    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))

    ## se
    expect_equal(diag(GS$cumhazard.se), test$cumhazard.se[,1])
    expect_equal(diag(GS$survival.se), test$survival.se[,1])
    
    ## iid
    GS.iid.diag <- do.call(rbind,lapply(1:NROW(dt),
                                        function(iN){GS$survival.iid[iN,iN,]}))
    expect_equal(GS.iid.diag, test$survival.iid[,1,])

    ## average.iid
    expect_equal(colMeans(GS.iid.diag), test2$survival.average.iid[,1])
    expect_equal(test$survival.average.iid, test2$survival.average.iid)

    ## average.iid with factor - diag=FALSE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(1:length(dt$time), nrow = NROW(dt), ncol = length(dt$time), byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = length(dt$time)))
    test3 <- predictCox(e.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)

    expect_equal(rowMultiply_cpp(GS$survival.average.iid, scale = 1:length(dt$time)), test3$survival.average.iid[[1]])
    expect_equal(t(apply(GS$survival.iid, 2:3, function(x){sum(x * (1:length(dt$time)))/length(x)})),
                 test3$survival.average.iid[[2]])

    ## average.iid with factor - diag=FALSE, time varying factor
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(rnorm(NROW(dt)*length(dt$time)), nrow = NROW(dt), ncol = length(dt$time)))
    test4 <- predictCox(e.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)
    expect_equal(do.call(rbind,lapply(1:NROW(dt), function(iObs){colMeans(GS$survival.iid[,,iObs] * attr(average.iid,"factor")[[1]])})),
                 test4$survival.average.iid[[1]])


    ## average.iid with factor - diag=TRUE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1, byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test5 <- predictCox(e.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = TRUE)

    
    expect_equal(5*test$survival.average.iid, test5$survival.average.iid[[1]])
    expect_equal(colMeans(colMultiply_cpp(GS.iid.diag, 1:length(dt$time))),
                 test5$survival.average.iid[[2]][,1])
    
})

test_that("[predictCox] diag strata", {
    eS.coxph <- coxph(Surv(time, event) ~ strata(X1) + X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCox(eS.coxph, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, average.iid = TRUE)
    test <- predictCox(eS.coxph, newdata = dt, times = dt$time,
                       se = FALSE, iid = TRUE, average.iid = TRUE, diag = TRUE)
    test2 <- predictCox(eS.coxph, newdata = dt, times = dt$time,
                       se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE)

    ## estimate
    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))

    ## iid
    GS.iid.diag <- do.call(rbind,lapply(1:NROW(dt),
                                        function(iN){GS$survival.iid[iN,iN,]}))
    expect_equal(GS.iid.diag, test$survival.iid[,1,])

    ## average.iid
    expect_equal(colMeans(GS.iid.diag), test2$survival.average.iid[,1])
    expect_equal(test$survival.average.iid, test2$survival.average.iid)

    ## average.iid with factor - diag=FALSE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(1:length(dt$time), nrow = NROW(dt), ncol = length(dt$time), byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = length(dt$time)))
    test3 <- predictCox(eS.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)

    expect_equal(rowMultiply_cpp(GS$survival.average.iid, scale = 1:length(dt$time)), test3$survival.average.iid[[1]])
    expect_equal(t(apply(GS$survival.iid, 2:3, function(x){sum(x * (1:length(dt$time)))/length(x)})),
                 test3$survival.average.iid[[2]])

    ## average.iid with factor - diag=FALSE, time varying factor
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(rnorm(NROW(dt)*length(dt$time)), nrow = NROW(dt), ncol = length(dt$time)))
    test4 <- predictCox(eS.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)
    expect_equal(do.call(rbind,lapply(1:NROW(dt), function(iObs){colMeans(GS$survival.iid[,,iObs] * attr(average.iid,"factor")[[1]])})),
                 test4$survival.average.iid[[1]])


    ## average.iid with factor - diag=TRUE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1, byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test5 <- predictCox(eS.coxph, newdata = dt, times = dt$time,
                        se = FALSE, iid = FALSE, average.iid = average.iid, diag = TRUE)

    
    expect_equal(5*test$survival.average.iid, test5$survival.average.iid[[1]])
    expect_equal(colMeans(colMultiply_cpp(GS.iid.diag, 1:length(dt$time))),
                 test5$survival.average.iid[[2]][,1])

})


## * [predictCox] Miscellaneous
## ** Confidence bands vs timereg
cat("[predictCox] Confidence band vs timereg \n")

## *** Data
set.seed(10)
dt <- sampleData(1e2, outcome = "survival")
newdata <- dt[1:10,]

dtStrata <- data.frame(time=c(4,3,1,1,2,2,3), 
                       status=c(1,1,1,0,1,1,0), 
                       x=c(0,2,1,1,1,0,0), 
                       sex=c(0,0,0,0,1,1,1)) 

## *** Model
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X2), data = dt, max.timepoint.sim=NULL)
e.coxph <- coxph(Surv(time, event) ~ X1 + X2, data = dt, x = TRUE, y = TRUE)

vec.times <- e.timereg$time.sim.resolution

## *** Compute quantile for confidence bands
resTimereg <- list()
for(i in 1:NROW(newdata)){ # i <- 1
    set.seed(10)
    resTimereg[[i]] <- predict.timereg(e.timereg,
                                       newdata = newdata[i,,drop=FALSE],
                                       times = vec.times,
                                       resample.iid = 1,
                                       n.sim = nsim.band)
}

## *** Tests
test_that("[predictCox] Quantile for the confidence band of the cumhazard", {

    predRR <- predictCox(e.coxph,
                         newdata = newdata,
                         times = vec.times,
                         se = TRUE,
                         iid = TRUE,
                         band = TRUE,
                         confint = FALSE,
                         type = "cumhazard")

    ## compatibility with timereg
    ref <- unlist(lapply(resTimereg,"[[", "unif.band"))

    set.seed(10)
    pred.band2 <- riskRegression:::confBandCox(iid = predRR$cumhazard.iid,
                                               se = predRR$cumhazard.se,
                                               n.sim = nsim.band, 
                                               conf.level = 0.95)

    expect_equal(pred.band2,ref)


    ## note confint is removing the first column since the standard error is 0
    set.seed(10)
    pred.band2.no0 <- riskRegression:::confBandCox(iid = predRR$cumhazard.iid[,-1,,drop=FALSE],
                                                   se = predRR$cumhazard.se[,-1],
                                                   n.sim = nsim.band, 
                                                   conf.level = 0.95)

    ## should not set transform to NA because at time 0 se=0 so the log-transform fails
    pred.confint <- confint(predRR, nsim.band = nsim.band, seed = 10,
                            cumhazard.transform = "none")
    expect_equal(pred.confint$cumhazard.quantileBand, pred.band2.no0)
    expect_equal(pred.confint$cumhazard.quantileBand, ref)

})

## *** Display
predRR <- predictCox(e.coxph,
                     newdata = newdata[1],
                     times = vec.times,
                     se = TRUE,
                     band = TRUE,
                     type = c("cumhazard","survival")
                     )


## plotRR <- autoplot(predRR, type = "survival", band = TRUE, ci = TRUE, plot = FALSE)
## dev.new()
## plotTR <- plot.predict.timereg(resTimereg[[1]])
## dev.new()
## plotRR$plot + coord_cartesian(ylim = c(0,1))
## graphics.off()

## *** With strata                                        
## Fit a stratified model 
eS.coxph <- coxph(Surv(time, status) ~ x + strata(sex), 
                  data = dtStrata, x = TRUE, y = TRUE) 

eS.pred  <- predictCox(eS.coxph, newdata = dtStrata, times = 1:4,
                       se = TRUE, iid = TRUE, band = TRUE)
eS.confint <- confint(eS.pred, seed = 10)
eS.confint$survival.quantileBand


## ** Dependence on data
cat("[predictCox] Dependence on data \n")
data(Melanoma, package = "riskRegression")

test_that("[predictCox] Dependence on data", {   
  Melanoma$entry <- -abs(rnorm(NROW(Melanoma), mean = 1, sd = 1))
  Melanoma2 <- Melanoma
  
  fit1 <- coxph(Surv(time, status>0)~strata(epicel)+age+strata(invasion)+sex+logthick, 
                data = Melanoma2, x = TRUE, y = TRUE)
  GS <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  
  Melanoma2 <- 7
  
  test <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  expect_equal(GS,test)
  
  ## with delayed entry
  Melanoma2 <- Melanoma
  
  fit1 <- coxph(Surv(entry ,time, status>0)~strata(epicel)+age+strata(invasion)+sex+logthick, 
                data = Melanoma2, x = TRUE, y = TRUE)
  GS <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  
  Melanoma2 <- 7
  
  test <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  expect_equal(GS,test)

})

## ** Store.iid argument
cat("[predictCox] Check same result store.iid=minimal vs. full \n")

## *** Data
set.seed(10)
d <- sampleData(50, outcome = "survival")
setkey(d,time)

## *** no strata
m.coxph <- coxph(Surv(time, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

## system.time(
##     res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata, store.iid = "minimal", se = TRUE, iid = FALSE)
## )
## system.time(
##     res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata, store.iid = "full", se = TRUE, iid = FALSE)
## )

test_that("[predictCox] store.iid = minimal vs. full - no strata", {
    res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", se = TRUE, iid = TRUE) 
    res2 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", average.iid = TRUE) 
    res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$cumhazard.se,res3$cumhazard.se)
    expect_equal(res1$survival.se,res3$survival.se)
    expect_equal(res1$cumhazard.iid,res3$cumhazard.iid)
    expect_equal(res1$survival.iid,res3$survival.iid)

    expect_equal(res2$cumhazard.average.iid, t(apply(res3$cumhazard.iid,2:3,mean)))
    expect_equal(res2$survival.average.iid, t(apply(res3$survival.iid,2:3,mean)))
})

## *** strata
m.coxph <- coxph(Surv(time, event) ~ strata(X1)+X6, data = d, y = TRUE, x = TRUE)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

test_that("[predictCox] store.iid = minimal vs. full - strata", {
    newdata <- rbind(d[1],d[1])
    res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", se = TRUE, iid = TRUE) 
    res2 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", average.iid = TRUE) 
    res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$cumhazard.se,res3$cumhazard.se)
    expect_equal(res1$survival.se,res3$survival.se)
    expect_equal(res1$cumhazard.iid,res3$cumhazard.iid)
    expect_equal(res1$survival.iid,res3$survival.iid)

    expect_equal(res2$cumhazard.average.iid, t(apply(res3$cumhazard.iid,2:3,mean)))
    expect_equal(res2$survival.average.iid, t(apply(res3$survival.iid,2:3,mean)))
})

## ** Weigthed cox
cat("[predictCox] Does not handle weights \n")
## *** Data
set.seed(10)
data(Melanoma, package = "riskRegression")
wdata <- runif(nrow(Melanoma), 0, 1)
times1 <- unique(Melanoma$time)

## *** Test
test_that("[predictCox] - weights",{

    fitW.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, weights = wdata, y = TRUE, x = TRUE)
    fitW.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE, weights = wdata)

    expect_error(resW <- predictCox(fitW.coxph, times = Melanoma$time, newdata = Melanoma))
    expect_error(resW <- predictCox(fitW.cph, times = Melanoma$time, newdata = Melanoma))
})

## ** Deal with negative timepoints
cat("[predictCox] dealing with negative timepoints or NA \n")
data(Melanoma, package = "riskRegression")

fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)

test_that("Deal with negative/NA time points",{
    expect_equal(unname(predictCox(fit.coxph, times = -1, newdata = Melanoma)$survival),
                 matrix(1,nrow = NROW(Melanoma), ncol = 1))

    expect_error(predictCox(fit.coxph, times = c(1,2,NA), newdata = Melanoma))
})
# }}}

## * [predictCox] Previous Bug
cat("[predictCox] Previous bug \n")
## ** Some coef are NA
dt <- sampleData(5e2, outcome = "survival")
e.coxph <- coxph(Surv(time, event) ~ X1+ X6 , data = dt, y = TRUE, x = TRUE)
e.coxph$coefficients[] <- as.numeric(NA)

test_that("Return error when coef contains NA", {
    expect_error(predictCox(e.coxph, newdata = dt, times = 1))
})

## ** average.iid

test_that("Cox - output of average.iid should not depend on other arguments", {
    set.seed(10)
    d <- sampleData(70,outcome="survival")
    d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]

    fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
                 data=d, ties="breslow", x = TRUE, y = TRUE)

    out1 <- predictCox(fit, newdata = d[1:5], times = 1:3, average.iid = TRUE)
    out2 <- predictCox(fit, newdata = d[1:5], times = 1:3, se = TRUE, average.iid = TRUE)

    expect_equal(out1$survival.average.iid,out2$survival.average.iid, tol = 1e-8)
})    


## ** Incorrect calculation of the standard error with atanh  (i.e. se/(1+b^2) instead of se/(1-b^2))
## from: Paul Blanche &lt;pabl@sund.ku.dk&gt;
## subject: suspected error in riskRegression
## date: Tue, 30 Jul 2019 11:42:14 +0200

test_that("Standard error after atanh transformation", {
    set.seed(10)
    x <- rnorm(1e2)
    y <- rnorm(1e2)

    rho <- cor.test(x,y)$estimate
    rho.se <- (1-rho^2)
    
    
    expect_equal(1, as.double(transformSE(estimate = rho, se = rho.se, type = "atanh")),
                 tol = 1e-5)
})



## ** se/iid should not depend on the ordering of the argument times

test_that("Cox - iid/se should not depend on other arguments", {
    set.seed(10)
    d <- sampleData(70,outcome="survival")
    d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]

    fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
                 data=d, ties="breslow", x = TRUE, y = TRUE)

    seqTau <- abs(rnorm(10))
    out1 <- predictCox(fit, newdata = d[1:5], times = seqTau,
                       se = TRUE, iid = TRUE, average.iid = TRUE)
    out2 <- predictCox(fit, newdata = d[1:5], times = sort(seqTau),
                       se = TRUE, iid = TRUE, average.iid = TRUE)

    out3 <- predictCox(fit, newdata = d[1:5], times = seqTau,
                       average.iid = TRUE)
    out4 <- predictCox(fit, newdata = d[1:5], times = sort(seqTau),
                       average.iid = TRUE)

    out5 <- predictCox(fit, newdata = d[1:5], times = seqTau,
                       se = TRUE, iid = TRUE, average.iid = TRUE, store.iid = "minimal")
    out6 <- predictCox(fit, newdata = d[1:5], times = sort(seqTau),
                       se = TRUE, iid = TRUE, average.iid = TRUE, store.iid = "minimal")

    expect_equal(out1$survival.iid[,order(seqTau),],out2$survival.iid)
    expect_equal(out1$survival.se[,order(seqTau)],out2$survival.se)
    expect_equal(out1$survival.average.iid[,order(seqTau)],out2$survival.average.iid)

    expect_equal(out2$survival.average.iid,out4$survival.average.iid)
    expect_equal(out3$survival.average.iid[,order(seqTau)],out4$survival.average.iid)

    expect_equal(out5$survival.iid[,order(seqTau),],out6$survival.iid)
    expect_equal(out5$survival.se[,order(seqTau)],out6$survival.se)
    expect_equal(out5$survival.average.iid[,order(seqTau)],out6$survival.average.iid)

    expect_equal(out2$survival.iid,out6$survival.iid)
    expect_equal(out2$survival.se,out6$survival.se)
    expect_equal(out2$survival.average.iid,out6$survival.average.iid)
})    


#----------------------------------------------------------------------
### test-predictCox.R ends here
