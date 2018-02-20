library(riskRegression)
library(testthat)
library(rms)
library(survival)

data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

## no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)

# {{{ 1- Check format of the output
cat("Format of the output \n")

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
# }}}

# {{{ 2- Check internal consistency
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

# {{{ 2- Check dependence on data
test_that("Dependence on data", {   
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
# }}}