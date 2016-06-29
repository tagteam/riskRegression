library(testthat)
library(riskRegression)
library(survival)
library(prodlim)
library(rms)
context("hazard/cumulative hazard/survival from Cox model")

#### data ####
data(Melanoma)

times1 <- unique(Melanoma$time[1:10])
times2 <- c(rnorm(10, mean = mean(Melanoma$time), sd = sd(Melanoma$time)),1e8)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

#### no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard - correct number of events",{ # does not work if strata because of identical event time in different strata
  lengthRes <- unlist(lapply(predictCox(fit.coxph, keep.times = TRUE), length))
  expect_equal(unname(lengthRes), rep(length(unique(fit.coxph$y[,"time"])), 4))
  lengthRes <- unlist(lapply(predictCox(fit.cph, keep.times = TRUE), length))
  expect_equal(unname(lengthRes), rep(length(unique(fit.cph$y[,"time"])), 4))
})


test_that("Refuse negative time points",{
  expect_error(predictCox(fit.coxph, times = -1, newdata = dataset1))
})

test_that("Prediction with Cox model - NA after last event",{
    test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
    prediction <- predictCox(fit.coxph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
    
    expect_equal(as.vector(is.na(prediction$hazard)), c(FALSE, FALSE, TRUE))
    expect_equal(as.vector(is.na(prediction$cumHazard)), c(FALSE, FALSE, TRUE))
    expect_equal(as.vector(is.na(prediction$survival)), c(FALSE, FALSE, TRUE))
  
})


#### strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard - order of the results",{
  expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), as.numeric(basehaz(fit.coxph)$strata))
  expect_equal(as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata)[predictCox(fit.cph, keep.strata = TRUE)$hazard>0], as.numeric(basehaz(fit.cph)$strata))
})

test_that("baseline hazard - match basehaz results",{
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumHazard, basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumHazard, basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
})

test_that("baseline hazard - maxtime",{
  # anywhere
  maxtime <- 3000
  predictTempo <- predictCox(fit.coxph, maxtime = maxtime, keep.times = TRUE)
  expect_equal(predictTempo$cumHazard, 
               basehaz(fit.coxph)$hazard[basehaz(fit.coxph)$time<maxtime], tolerance = 1e-8)
  
  # at an event time
  maxtime <- Melanoma$time[100]
  predictTempo <- predictCox(fit.coxph, maxtime = maxtime, keep.times = TRUE)
  expect_equal(predictTempo$cumHazard[predictTempo$time!=maxtime], 
               basehaz(fit.coxph)$hazard[basehaz(fit.coxph)$time<maxtime], tolerance = 1e-8)
})

test_that("Prediction with Cox model - export of strata and times",{
  predictTempo <- predictCox(fit.coxph)
  expect_equal(length(predictTempo$strata)>0, FALSE) 
  expect_equal(length(predictTempo$times)>0, FALSE)
  predictTempo <- predictCox(fit.coxph, keep.strata = TRUE) # as.data.table(predictCox(fit.coxph, keep.strata = TRUE))
  expect_equal(length(predictTempo$strata)>0, TRUE) 
  expect_equal(length(predictTempo$times)>0, FALSE)
  predictTempo <- predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE)
  expect_equal(length(predictTempo$strata)>0, TRUE)
  expect_equal(length(predictTempo$times)>0, TRUE)
  
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1)
  expect_equal(length(predictTempo$strata)>0, FALSE) 
  expect_equal(length(predictTempo$times)>0, FALSE)
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1, keep.strata = TRUE)
  expect_equal(length(predictTempo$strata)>0, TRUE)
  expect_equal(length(predictTempo$times)>0, FALSE)
  predictTempo <- predictCox(fit.coxph, times = sort(times2), newdata = dataset1, keep.strata = TRUE, keep.times = TRUE)
  expect_equal(length(predictTempo$strata)>0, TRUE)
  expect_equal(length(predictTempo$times)>0, TRUE)
})

test_that("Prediction with Cox model - consistency of hazard/cumHazard/survival",{
  predictTempo <- predictCox(fit.coxph, times = times1, newdata = dataset1)
  expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumHazard,1,diff)), tolerance = 1e-8)
  expect_equal(predictTempo$survival, exp(-predictTempo$cumHazard), tolerance = 1e-8)
})

#### test correct prediction after the last event 
# take one observation from each strata
data.test <- data.table(Melanoma)[, .SD[1], by = c("invasion", "ici")]
setkeyv(data.test, c("invasion","ici"))

# identify the last event time for each strata
epsilon <- min(diff(unique(fit.coxph$y[,"time"])))/10
baseHazStrata <- as.data.table(predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE))
dt.times <- baseHazStrata[, .(beforeLastTime = times[.N]-epsilon, LastTime = times[.N], afterLastTime = times[.N]+epsilon), by = strata]

test_that("Prediction with Cox model - NA after last event",{
  for(Ttempo in 1:nrow(dt.times)){
  test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
  
  prediction <- predictCox(fit.coxph, times = test.times, newdata = data.test)
  expect_equal(is.na(prediction$hazard[Ttempo,]), c(FALSE, FALSE, TRUE))
  expect_equal(is.na(prediction$cumHazard[Ttempo,]), c(FALSE, FALSE, TRUE))
  expect_equal(is.na(prediction$survival[Ttempo,]), c(FALSE, FALSE, TRUE))
  }
  # pec:::predictSurvProb(fit.coxph, times = test.times, newdata = data.test)
})



