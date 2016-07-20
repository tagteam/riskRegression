### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:35) 
## Version: 
## last-updated: Jun  6 2016 (16:52) 
##           By: Thomas Alexander Gerds
##     Update #: 16
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(riskRegression)
library(testthat)
library(rms)
library(survival)
context("Risk prediction")

#### 1- Comparison between results from riskRegression and pec package ####

#### CSC  model
## predictRisk (predict.CauseSpecificCox function) vs. predictEventProb
test_that("Extract risks from coxph and cph and CSC objects",{
    set.seed(10)
    n <- 300
    df <- SimCompRisk(n)
    dn <- SimCompRisk(17)
    df$time <- round(df$time,2)
    F1 <- CSC(Hist(time,event) ~ X1+X2,data = df,x = TRUE, y = TRUE )
    ttt <- sort(unique(df$time))
    ## cause 1
    a <- pec::predictEventProb(F1,times=ttt,newdata=dn,cause=1)
    b <- riskRegression::predictRisk(F1,times=ttt,newdata=dn,cause=1)
    attributes(a) <- attributes(b)
    expect_equal(a,b)
    ## cause 2
    a <- pec::predictEventProb(F1,times=ttt,newdata=dn,cause=2)
    b <- riskRegression::predictRisk(F1,times=ttt,newdata=dn,cause=2)
    attributes(a) <- attributes(b)
    expect_equal(a,b)
    ## Cox object 1
    a <- pec::predictSurvProb(F1$models[[1]],times=ttt,newdata=dn)
    b <- 1-riskRegression::predictRisk(F1$model[[1]],times=ttt,newdata=dn)
    attributes(a) <- attributes(b)
    expect_equal(a,b)
    ## Cox object 2
    a <- pec::predictSurvProb(F1$models[[2]],times=ttt,newdata=dn)
    b <- 1-riskRegression::predictRisk(F1$model[[2]],times=ttt,newdata=dn)
    attributes(a) <- attributes(b)
    expect_equal(a,b)
})

## fitter and method to handle ties - no strata
set.seed(10)
n <- 300
df <- SimCompRisk(n)
dn <- SimCompRisk(17)
seqTime <- c(unique(sort(df$time)), max(df$time) + 1)[1:10] # to be removed in future versions
for(model in c("coxph","cph")){
    for(method.ties in c("breslow","efron")){
        CSC.NS <- CSC(Hist(time,event) ~ X1*X2,data = df, method = method.ties, fitter = model)
        for(cause in 1:2){ 
            Event0.NS <- predictRisk(CSC.NS, newdata = dn, times = seqTime, cause = cause)
            EventTest.NS <- pec::predictEventProb(CSC.NS, newdata = dn, times = seqTime, cause = cause)
            test_that(paste0("predictRisk vs predictEventProb - no strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
                expect_equal(as.numeric(Event0.NS),as.numeric(unlist(EventTest.NS)), tolerance = 1e-8)
            })
        }
    }
}

## fitter and method to handle ties - with strata
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[1:10] # to be removed in future versions
for(model in c("coxph","cph")){
    for(method.ties in c("breslow","efron")){
        if(model == "coxph"){
            CSC.S <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties, fitter = "coxph")
        }else if(model == "cph"){
            CSC.S <- CSC(Hist(time,event) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph" )
        }
        for(cause in 1:2){ 
            Event0.S <- predictRisk(CSC.S, newdata = df.S, times = seqTime, cause = cause)
            EventTest.S <- pec::predictEventProb(CSC.S, newdata = df.S, times = seqTime, cause = cause)
            test_that(paste0("predictRisk vs predictEventProb - strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
                expect_equal(as.numeric(Event0.S),as.numeric(unlist(EventTest.S)), tolerance = 1e-8)
            })
        }
    }
}

#### Cox model: predictRisk (predictCox function) vs. predictSurvProb
## no strata
set.seed(10)
n <- 100
df <- SimSurv(n)
df$time <- round(df$time,1)
## ttt <- c(seq(0, max(df$time), length.out = 10), max(df$time) + 1)
ttt <- c(seq(0, max(df$time)/2, length.out = 10))
for(model in c("coxph","cph")){
    for(method.ties in c("breslow","efron")){
        if(model == "coxph"){
            coxph <- coxph(formula = Surv(time,status) ~ X1 +  X2, data = df, ties = method.ties)
            Surv0 <- predictRisk(coxph, newdata = df, times = ttt)
        }else{
            coxph <- cph(formula = Surv(time,status) ~ X1 + X2, data = df, ties = method.ties, surv = TRUE, x = TRUE, y = TRUE)
            Surv0 <- predictRisk(coxph, newdata = df, times = ttt)
        }
        SurvProb <- pec::predictSurvProb(coxph, ttt, newdata = df)
        test_that(paste0("predictRisk vs predictSurvProb - without strata, method.ties/method.baseHaz/model = ",method.ties,"/",model,""),{
            expect_equal(object = Surv0, expected = 1-SurvProb, tolerance=1e-8)
        })
    }
}

## strata
set.seed(10)
n <- 100
df <- SimSurv(n)
df$time <- round(df$time,1)
df$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
## ttt <- c(seq(0, max(df$time), length.out = 10), max(df$time) + 1)
ttt <- c(seq(0, max(df$time), length.out = 10))
for(model in c("coxph","cph")){
  ## print(model)
  for(method.ties in c("breslow","efron")){
    ## print(method.ties)
    if(model == "coxph"){
      fit <- coxph(formula = Surv(time,status) ~ strata(X1) + strata(X3) + X2, data = df, method = method.ties)
      Surv0 <- predictRisk(fit, newdata = df, times = ttt)
    }else{
      fit <- cph(formula = Surv(time,status) ~ strat(X1) + strat(X3) + X2, data = df, method = method.ties, surv = TRUE, x = TRUE, y = TRUE)
      Surv0 <- predictRisk(fit, newdata = df, times = ttt)
    }
    SurvProb <- pec::predictSurvProb(fit, ttt, newdata = df)
    test_that(paste0("predictRisk vs predictSurvProb - without strata, method.ties/model = ",method.ties,"/",model,""),{
      a <- Surv0
      b <- 1-SurvProb
      attributes(a) <- attributes(b)
      expect_equal(object = a, expected = b, tolerance=1e-6)
    })
  }
}

#### 2- Benefit of using maxtime ####
# limits the computation of the baseline hazard to time<=maxtime

set.seed(10)
n <- 1e5
df <- SimCompRisk(n)
dn <- SimCompRisk(17)
F1 <- CSC(Hist(time,event) ~ X1+X2,data = df)

timeBase <- system.time(
  resBase <- predictCox(F1$models[[1]], newdata=dn, times = c(0,0.1), maxtime = Inf)
)
timeMaxTime <- system.time(
  resMaxTime <- predictCox(F1$models[[1]], newdata=dn, times = c(0,0.1))
)
print( (timeBase-timeMaxTime) )
test_that("baseline hazard (no strata) - maxtime",{ 
  expect_equal(resBase, resMaxTime)
})


F1 <- CSC(Hist(time,event) ~ strata(X1)+X2, data = df)

timeBase <- system.time(
  resBase <- predictCox(F1$models[[1]], newdata=dn, times = c(0,0.1), maxtime = Inf)
)
timeMaxTime <- system.time(
  resMaxTime <- predictCox(F1$models[[1]], newdata=dn, times = c(0,0.1))
)
print( (timeBase-timeMaxTime) )

test_that("baseline hazard (strata) - maxtime",{ 
  expect_equal(resBase, resMaxTime)
})

#### 
n <- 3
set.seed(3)
dn <- SimCompRisk(n)
dn$time <- round(dn$time,2)
dn$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
CSC.h3 <- CSC(Hist(time,event) ~ X1 + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph" )
CSC.h1 <- CSC(Hist(time,event) ~ strat(X1) + X3 + X2, data = df.S, ties = method.ties, fitter = "cph" )
CSC.h <- CSC(Hist(time,event) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph" )
CSC.s <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties, fitter = "coxph" )
predictRisk(CSC.h1, newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.h3, newdata = dn, times = c(5,10,15,20), cause = cause)

CSC.h0 <- CSC(Hist(time,event) ~ X1 + X3 + X2, data = df.S, ties = method.ties, fitter = "cph" )
predictRisk(CSC.h0, newdata = dn, times = c(5,10,15,20), cause = cause)

predictRisk(CSC.h1, newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s, newdata = dn, times = c(5,10,15,20), cause = cause)


df.S[df.S$time==6.55,c("time","event")]
predictCox(CSC.h$models[[1]],newdata = dn[1,],times=c(2.29,6.55),type="hazard")
predictCox(CSC.h$models[[2]],newdata = dn[1,],times=6.55,type="hazard")

predictCox(CSC.h$models[["Cause 1"]],newdata = dn[1,],times=CSC.h$eventTimes,type="hazard")$hazard

predictRisk(CSC.h, newdata = dn[1,], times = c(2), cause = "1")
predictRisk(CSC.h, newdata = dn, times = c(1,2,3.24,3.25,3.26,5,10,15,20), cause = cause)
predictRisk(CSC.s, newdata = dn, times = c(1,5,10,15,20), cause = cause)

predictCox(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20))

predictRisk(CSC.h$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)

predictRisk(CSC.h$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)


#### 3- Check prediction after and before the last event ####
data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
fit.CSC <- CSC(Hist(time,event) ~ thick + strat(invasion) + strat(ici), data = Melanoma, fitter = "cph")

# take one observation from each strata
data.test <- data.table(Melanoma)[, .SD[1], by = c("invasion", "ici")]
setkeyv(data.test, c("invasion","ici"))

# identify the last event time for each strata
epsilon <- min(diff(unique(fit.coxph$y[,"time"])))/10
baseHazStrata <- as.data.table(predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE))
dt.times <- baseHazStrata[, .(beforeLastTime = times[.N]-epsilon, LastTime = times[.N], afterLastTime = times[.N]+epsilon), by = strata]

#### predictCox
test_that("Prediction with Cox model (strata) - NA after last event",{
  for(Ttempo in 1:nrow(dt.times)){
    test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
    
    prediction <- predictCox(fit.coxph, times = test.times, newdata = data.test)
    expect_equal(is.na(prediction$hazard[Ttempo,]), c(FALSE, FALSE, TRUE))
    expect_equal(is.na(prediction$cumHazard[Ttempo,]), c(FALSE, FALSE, TRUE))
    expect_equal(is.na(prediction$survival[Ttempo,]), c(FALSE, FALSE, TRUE))
    
    prediction <- predictCox(fit.cph, times = test.times, newdata = data.test)
    expect_equal(is.na(prediction$hazard[Ttempo,]), c(FALSE, FALSE, TRUE))
    expect_equal(is.na(prediction$cumHazard[Ttempo,]), c(FALSE, FALSE, TRUE))
    expect_equal(is.na(prediction$survival[Ttempo,]), c(FALSE, FALSE, TRUE))
  }
  
})

test_that("Prediction with CSC (strata) - NA after last event",{
  for(Ttempo in 1:nrow(dt.times)){
    test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
    
    prediction <- predict(fit.CSC, times = test.times, newdata = data.test, cause = "death.malignant.melanoma")
    expect_equal(is.na(prediction[Ttempo,]), c(FALSE, FALSE, TRUE))
    }
  
})

#### 4- Others ####
data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

## no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard - correct number of events",{ 
  lengthRes <- unlist(lapply(predictCox(fit.coxph, keep.times = TRUE), length))
  expect_equal(unname(lengthRes), rep(length(unique(fit.coxph$y[,"time"])), 4))
  lengthRes <- unlist(lapply(predictCox(fit.cph, keep.times = TRUE), length))
  expect_equal(unname(lengthRes), rep(length(unique(fit.cph$y[,"time"])), 4))
})

# test_that("Refuse negative time points",{
#   expect_error(predictCox(fit.coxph, times = -1, newdata = dataset1))
# })

test_that("Prediction with Cox model - NA after last event",{
  test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
  prediction <- predictCox(fit.coxph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  
  expect_equal(as.vector(is.na(prediction$hazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$cumHazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$survival)), c(FALSE, FALSE, TRUE))
})

test_that("Prediction with Cox model - sorted vs. unsorted times",{
  newOrder <- sample.int(length(times2),length(times2),replace = FALSE)
  predictionUNS <- predictCox(fit.coxph, times = times2[newOrder], newdata = Melanoma)
  predictionS <- predictCox(fit.coxph, times = times2, newdata = Melanoma)
  
  expect_equal(predictionS, lapply(predictionUNS, function(x){x[,order(newOrder)]}))
})

# test_that("Prediction with Cox model - NA in times",{
#   predictionS <- predictCox(fit.coxph, times = c(times2,NA), newdata = Melanoma)
# })

## strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard (strata) - order of the results",{
  expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), as.numeric(basehaz(fit.coxph)$strata))
  expect_equal(as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata)[predictCox(fit.cph, keep.strata = TRUE)$hazard>0], as.numeric(basehaz(fit.cph)$strata))
})

test_that("baseline hazard (strata) - match basehaz results",{
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumHazard, basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumHazard, basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
})

test_that("baseline hazard (strata) - correct number of events",{ 
  strata <- interaction(Melanoma$invasion, Melanoma$ici)
  timePerStrata <- tapply(fit.coxph$y[,"time"],strata, function(x){length(unique(x))})
  
  lengthRes <- unlist(lapply(predictCox(fit.coxph, keep.times = TRUE, keep.strata = TRUE), length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
  lengthRes <- unlist(lapply(predictCox(fit.cph, keep.times = TRUE, keep.strata = TRUE), length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
})

test_that("baseline hazard (strata) - maxtime",{
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

test_that("Prediction with Cox model (strata) - export of strata and times",{
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

test_that("Prediction with Cox model (strata) - consistency of hazard/cumHazard/survival",{
  predictTempo <- predictCox(fit.coxph, times = times1, newdata = dataset1)
  expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumHazard,1,diff)), tolerance = 1e-8)
  expect_equal(predictTempo$survival, exp(-predictTempo$cumHazard), tolerance = 1e-8)
})


test_that("Prediction with Cox model (strata) - sorted vs. unsorted times",{
  newOrder <- sample.int(length(times2),length(times2),replace = FALSE)
  predictionUNS <- predictCox(fit.coxph, times = times2[newOrder], newdata = Melanoma)
  predictionS <- predictCox(fit.coxph, times = times2, newdata = Melanoma)
  
  expect_equal(predictionS, lapply(predictionUNS, function(x){x[,order(newOrder)]}))
})

test_that("Prediction with Cox model (strata) - incorrect strata",{
  dataset1$invasion <- "5616"
  expect_error(predictCox(fit.coxph, times = times1, newdata = dataset1))
})

test_that("Prediction with Cox model (strata) - no event before prediction time",{
  data(lung, package = "survival")
  fitS <- coxph(Surv(time, status) ~ age + strata(ph.ecog) + inst, lung[1:20,])
  predictCox(fitS, newdata = lung[1:10,], times = 10)
})

#----------------------------------------------------------------------
### predictRisk.R ends here
