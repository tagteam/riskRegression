### test-predict-SEconfint.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: sep  8 2019 (16:35) 
##           By: Brice Ozenne
##     Update #: 197
#----------------------------------------------------------------------
## 
### Commentary: 
## Compare the survival and its standard error obtained with iidCox and timereg:
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Setting
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(timereg)

context("[predictCox] Computation of iid,SE,CI,CB, comparison to timereg")
nsim.band <- 500

## * predictions: no strata
cat("[predictCox] Predictions (no strata, continuous) \n")

## ** Data
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]
dt[,X1:=as.numeric(as.character(X1))]
dt[,X2:=as.numeric(as.character(X2))]
dt[ , X16 := X1*X6]
dt[ , Xcat2 := as.factor(paste0(X1,X2))]

## sorted dataset
dt.sort <- copy(dt)
setkeyv(dt.sort,c("time")) 

## ** Model
e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)
e.cph <- cph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)
e.coxph_sort <- coxph(Surv(time, event) ~ X1*X6, data = dt.sort, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X6) + prop(X1*X6),
                       data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

## ** consistency between hazard/cumhazard/survival

test_that("[predictCox] - consistency of hazard/cumhazard/survival",{
  predRR <- predictCox(e.coxph, type = c("hazard","cumhazard","survival"), times = sort(dt$time), newdata = dt)
  expect_equal(predRR$hazard[,-1], t(apply(predRR$cumhazard,1,diff)), tolerance = 1e-8)
  expect_equal(predRR$survival, exp(-predRR$cumhazard), tolerance = 1e-8)
})

## ** 1 fixed time
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

## ** after last event
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


## ** before first event
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

## ** sorted vs. unsorted times

test_that("[predictCox] - sorted vs. unsorted times (no strata)",{
    vec.time <- dt$time
    index.sort <- order(vec.time)
    vec.time_sort <- dt$time[index.sort]
    
    predRR1 <- predictCox(e.coxph, newdata = dt, times = vec.time, se = TRUE)
    predRR2 <- predictCox(e.coxph, newdata = dt, times = vec.time_sort, se = TRUE)

    expect_equal(predRR1$survival[,index.sort],
                 predRR2$survival)

    expect_equal(predRR1$time[index.sort],
                 predRR2$time)

})

## ** iid.average
test_that("[predictCox] - fast iid average (no strata)",{
    ## simple average
    predRR.av <- predictCox(e.coxph, times = dt$time[1:5], average.iid = TRUE, newdata = dt,
                            type = c("cumhazard","survival"))
    predRR.GS <- predictCox(e.coxph, times = dt$time[1:5], iid = TRUE, newdata = dt,
                            type = c("cumhazard","survival"))

      
    expect_equal(t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$cumhazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$survival.iid, MARGIN = 2:3,mean)),
                 predRR.av$survival.average.iid, tolerance = 1e-8)

    ## weighted average
    fT <- TRUE
    attr(fT, "factor") <- list(matrix(1, nrow = NROW(dt), ncol = 5),
                               matrix(1:NROW(dt), nrow = NROW(dt), ncol = 5)
                               )
    
    predRR.av2 <- predictCox(e.coxph, times = sort(dt$time[1:5]), average.iid = fT, newdata = dt,
                             type = c("cumhazard","survival"))
    GS <- t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3, function(iCol){
        mean(iCol * attr(fT, "factor")[[2]][,1])
    }))

    expect_equal(predRR.av$cumhazard.average.iid[,order(dt$time[1:5])],
                 predRR.av2$cumhazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(GS[,order(dt$time[1:5])],
                 predRR.av2$cumhazard.average.iid[[2]],
                 tolerance = 1e-8)
})

## * predictions: no strata with a categorical variable
cat("[predictCox] Predictions (no strata, categorical) \n")
## ** Data
## see no strata
## ** Model
e.coxph <- coxph(Surv(time, event) ~ Xcat2 + X6 , data = dt, y = TRUE, x = TRUE)
e.cph <- cph(Surv(time, event) ~ Xcat2 + X6 , data = dt, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(Xcat2) + prop(X6),
                       data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

## ** test


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

## * predictions: strata
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
                            type = c("cumhazard","survival"))
    predRR.GS <- predictCox(eS.coxph, times = seqTime, iid = TRUE, newdata = dtStrata,
                            type = c("cumhazard","survival"))

      
    expect_equal(t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3,mean)),
                 predRR.av$cumhazard.average.iid, tolerance = 1e-8)

    expect_equal(t(apply(predRR.GS$survival.iid, MARGIN = 2:3,mean)),
                 predRR.av$survival.average.iid, tolerance = 1e-8)

    
    ## only one strata
    index.strata1 <- dtStrata[,.I[strata == 1]]
    predRR.avStrata1 <- predictCox(eS.coxph, times = seqTime, average.iid = TRUE, newdata = dtStrata[index.strata1],
                                   type = c("cumhazard","survival"))

    expect_equal(t(apply(predRR.GS$cumhazard.iid[index.strata1,,], MARGIN = 2:3,mean)),
                 predRR.avStrata1$cumhazard.average.iid, tolerance = 1e-8)
    expect_equal(t(apply(predRR.GS$survival.iid[index.strata1,,], MARGIN = 2:3,mean)),
                 predRR.avStrata1$survival.average.iid, tolerance = 1e-8)

    ## weighted average
    fT <- TRUE
    attr(fT, "factor") <- list(matrix(1, nrow = NROW(dtStrata), ncol = length(seqTime)),
                               matrix(1:NROW(dtStrata), nrow = NROW(dtStrata), ncol = length(seqTime))
                               )
    
    predRR.av2 <- predictCox(eS.coxph, times = sort(seqTime), average.iid = fT, newdata = dtStrata,
                             type = c("cumhazard","survival"))
    GS <- t(apply(predRR.GS$cumhazard.iid, MARGIN = 2:3, function(iCol){
        mean(iCol * attr(fT, "factor")[[2]][,1])
    }))

    expect_equal(predRR.av$cumhazard.average.iid[,order(seqTime)],
                 predRR.av2$cumhazard.average.iid[[1]],
                 tol = 1e-8)
    expect_equal(GS[,order(seqTime)],
                 predRR.av2$cumhazard.average.iid[[2]],
                 tolerance = 1e-8)
})

## * predictions: SE/CI check against manual computation
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


## * Dependence on data
cat("[predictCox] Dependence on data \n")
data(Melanoma)

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

## * Store.iid argument
cat("[predictCox] Check same result store.iid=minimal vs. full \n")

## ** Data
set.seed(10)
d <- sampleData(50, outcome = "survival")
setkey(d,time)

## ** no strata
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

## ** strata
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

# }}}

## * Weigthed cox
cat("[predictCox] does not handle weights \n")
## ** Data
set.seed(10)
data(Melanoma)
wdata <- runif(nrow(Melanoma), 0, 1)
times1 <- unique(Melanoma$time)

## ** Test
test_that("[predictCox] - weights",{

fitW.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, weights = wdata, y = TRUE, x = TRUE)
fitW.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE, weights = wdata)

expect_error(resW <- predictCox(fitW.coxph, times = Melanoma$time, newdata = Melanoma))
expect_error(resW <- predictCox(fitW.cph, times = Melanoma$time, newdata = Melanoma))
})
# }}}

## * Time varying variables
# NOTE: I am not sure what timereg is exactly doing in this case

## d <- sampleData(1e2, outcome = "survival")
## d$start <- runif(NROW(d),min=0,max=(d$eventtime-0.1) )

## test_that("predicted survival with time varying variables",{
    ## fit.coxph <- coxph(Surv(start, eventtime, event) ~ X1, data = d, x = TRUE, y = TRUE)
    ## fit.timereg <- cox.aalen(Surv(start, eventtime, event) ~ prop(X1), data = d)
    ## predTimes <- sort(unique(d$eventtime))
    ## M1 <- predictCox(fit.coxph, newdata = d, time = predTimes)$survival
    ## M2 <- predict(fit.timereg, newdata = d, time = predTimes)$S0
    ## expect_equal(M1,M2)
## })
# }}}

## * Diag argument
cat("[predictCox] diag argument \n")
set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]

test_that("[predictCox] diag no strata", {
    e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCox(e.coxph, newdata = dt, times = dt$time, se = FALSE, iid = TRUE)
    test <- predictCox(e.coxph, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, diag = TRUE)
    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))

    GS <- predictCox(e.coxph, newdata = dt, times = dt$time, se = FALSE, iid = TRUE)
    test <- predictCox(e.coxph, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, diag = TRUE)

    GS.iid.diag <- do.call(rbind,lapply(1:NROW(dt),
                                        function(iN){GS$cumhazard.iid[iN,iN,]}))
    expect_equal(GS.iid.diag, test$cumhazard.iid[,1,])
})

test_that("[predictCox] diag strata", {
    eS.coxph <- coxph(Surv(time, event) ~ strata(X1) + X6, data = dt, y = TRUE, x = TRUE)

    GS <- predictCox(eS.coxph, newdata = dt, times = dt$time, se = FALSE)
    test <- predictCox(eS.coxph, newdata = dt, times = dt$time, se = FALSE, diag = TRUE)

    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$cumhazard), as.double(test$cumhazard))
    expect_equal(diag(GS$survival), as.double(test$survival))
})

## * Previous Bug
cat("[predictCox] Previous bug \n")
## ** Some coef are NA
dt <- sampleData(5e2, outcome = "survival")
e.coxph <- coxph(Surv(time, event) ~ X1+ X6 , data = dt, y = TRUE, x = TRUE)
e.coxph$coefficients[] <- as.numeric(NA)

test_that("Return error when coef contains NA", {
    expect_error(predictCox(e.coxph, newdata = dt, times = 1))
})

## ** 0 average.iid

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

test_that("CSC - output of average.iid should not depend on other arguments", {
    set.seed(10)
    d <- sampleData(70,outcome="competing.risks")
    d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]

    fit <- CSC(Hist(time,event)~X1 + strata(X2) + X6,
               data=d)

    out1 <- predict(fit, newdata = d[1:5], times = 1:3, average.iid = TRUE, cause = 1)
    out2 <- predict(fit, newdata = d[1:5], times = 1:3, se = TRUE, average.iid = TRUE, cause = 1)

    test_that("output of average.iid should not depend on other arguments", {
        expect_equal(out1$survival.average.iid,out2$survival.average.iid, tol = 1e-8)
    })    
})

## ** (Previously) incorrect calculation of the standard error with atanh  (i.e. se/(1+b^2) instead of se/(1-b^2))
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

    out1 <- predictCox(fit, newdata = d[1:5], times = c(0,3), se = TRUE)
    out2 <- predictCox(fit, newdata = d[1:5], times = c(3,0), se = TRUE)

    expect_equal(out1$survival.se,out2$survival.se[,2:1], tol = 1e-8)

    ## d <- sampleData(70)
    ## d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]
    ## fit <- CSC(Hist(time,event) ~ X1 + X2, data = d)

    ## out1 <- predict(fit, newdata = d[1:5], times = c(0,3), se = TRUE, cause = 1)
    ## out2 <- predict(fit, newdata = d[1:5], times = c(3,0), se = TRUE, cause = 1)
})    



## ** ???
f1 <- coxph(Surv(time,status==1) ~ age+logthick+epicel+strata(sex),
            data=Melanoma, x=TRUE,y=TRUE)
res <- predictCox(f1,newdata=Melanoma[c(17,101,123),],
                  times=c(7,3,5)*365.25)
}

##----------------------------------------------------------------------
### test-predictCox-SEconfint.R ends here



## * Timing
if(FALSE){

    ## library(microbenchmark)
    ## library(survival)
    ## library(rms)
    ## library(mets)
    ## library(profvis)
    
    ## set.seed(10)
    ## d <- sampleData(1e3)
    ## m.coxph <- coxph(Surv(time, event>0) ~ X1 + X2 + X3, data = d, x = TRUE, y = TRUE)
    ## m.cph <- cph(Surv(time, event>0) ~ X1 + X2 + X3, data = d, x = TRUE, y = TRUE)
    ## m.phreg <- phreg(Surv(time, event>0) ~ X1 + X2 + X3, data = d)


    ## ## baseline hazard
    ## runBenchmark <- function(n) {
    ##     microbenchmark(times = 20,  
    ##                    coxph = {predictCox(m.coxph)},
    ##                    cph = {predictCox(m.cph)},
    ##                    phreg = {predictCox(m.phreg)},
    ##                    basehaz = {basehaz(m.coxph)}
    ##                    )
    ## }
    ## res <- runBenchmark()

    ## profvis(
    ##     predictCox(m.coxph)
    ## )

    ## ## predictions
    ## runBenchmark <- function(n) {
    ##     microbenchmark(times = 20,  
    ##                    coxph = {predictCox(m.coxph, newdata = d, times = 1:5)},
    ##                    cph = {predictCox(m.cph, newdata = d, times = 1:5)},
    ##                    phreg = {predictCox(m.phreg, newdata = d, times = 1:5)}
    ##                    )
    ## }
    ## res <- runBenchmark()
    ## res

    ## profvis(
    ##     predictCox(m.coxph, newdata = d, times = 1:5)
    ## )

    ## ## predictions with SE
    ## runBenchmark <- function(n) {
    ##     microbenchmark(times = 20,  
    ##                    coxph = {predictCox(m.coxph, newdata = d, times = 1:5, se = TRUE)},
    ##                    cph = {predictCox(m.cph, newdata = d, times = 1:5, se = TRUE)},
    ##                    phreg = {predictCox(m.phreg, newdata = d, times = 1:5, se = TRUE)}
    ##                    )
    ## }
    ## res <- runBenchmark()
    ## res

    ## profvis(
    ##     predictCox(m.coxph, newdata = d, times = 1:5, se = TRUE)
    ## )
}
