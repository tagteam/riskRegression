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
cat("riskRegression vs pec \n")

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
            expect_equal(object = as.double(Surv0), expected = as.double(1-SurvProb), tolerance=1e-8)
        })
        # as.double because predictRisk do not set colnames while predictSurvProb do set times as colnames
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


#### 2- [predictCox,CSC] Check prediction after and before the last event ####
cat("prediction before and after the last event \n")

data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

#### no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma, fitter = "cph")

test_that("Prediction with Cox model - NA after last event",{
  test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
  
  prediction <- predictCox(fit.coxph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.vector(is.na(prediction$hazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$cumHazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$survival)), c(FALSE, FALSE, TRUE))
  
  prediction <- predictCox(fit.cph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.vector(is.na(prediction$hazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$cumHazard)), c(FALSE, FALSE, TRUE))
  expect_equal(as.vector(is.na(prediction$survival)), c(FALSE, FALSE, TRUE))
})

test_that("Prediction with Cox model - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predictCox(fit.coxph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.double(prediction$hazard), 0)
  expect_equal(as.double(prediction$cumHazard), 0)
  expect_equal(as.double(prediction$survival), 1)
  
  prediction <- predictCox(fit.cph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.double(prediction$hazard), 0)
  expect_equal(as.double(prediction$cumHazard), 0)
  expect_equal(as.double(prediction$survival), 1)
})

test_that("Prediction with CSC - NA after last event",{
  test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.vector(is.na(prediction)), c(FALSE, FALSE, TRUE))
})


test_that("Prediction with CSC - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction), 0)
})

test_that("Prediction - last event censored",{
  nData <- length(Melanoma$event)
  
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma)
  fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
  fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma, fitter = "cph")

  vec.times <- sort(c(Melanoma$time, Melanoma$time+1/2)) 
  p1 <- predictCox(fit.coxph, times = vec.times, newdata = Melanoma)
  p2 <- predictCox(fit.cph, times = vec.times, newdata = Melanoma)
  p3 <- pec::predictSurvProb(fit.coxph, times = vec.times, newdata = Melanoma) # predictSurvProb automatically sort the results
  
  expect_equal(p1,p2, tolerance = 1e-4)
  expect_equal(p1$survival,p3)
  
  survLast <- p1$survival[,vec.times==Melanoma$time[nData]]
  survLastM1 <- p1$survival[,vec.times==Melanoma$time[nData-1]]
  expect_equal(unname(survLast-survLastM1==0), rep(TRUE, nData)) # check that survival decrease at the last time
})

test_that("Prediction - last event is a death",{
  nData <- length(Melanoma$event)
  Melanoma2 <- Melanoma
  Melanoma2$status[nData] <- 1
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma2)
  fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma2, y = TRUE, x = TRUE)
  fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma2, fitter = "cph")
  
  vec.times <-  sort(c(Melanoma$time, Melanoma$time+1/2)) 
  p1 <- predictCox(fit.coxph, times = vec.times, newdata = Melanoma)
  p2 <- predictCox(fit.cph, times = vec.times, newdata = Melanoma)
  p3 <- pec::predictSurvProb(fit.coxph, times = vec.times, newdata = Melanoma) # predictSurvProb automatically sort the results
  
  expect_equal(p1,p2, tolerance = 1e-4)
  expect_equal(p1$survival,p3)
  
  survLast <- p1$survival[,vec.times==Melanoma2$time[nData]]
  survLastM1 <- p1$survival[,vec.times==Melanoma2$time[nData-1]]
  expect_equal(unname(survLast-survLastM1<0), rep(TRUE, nData)) # check that survival decrease at the last time
})

#### strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion) + strat(ici), data = Melanoma, fitter = "cph")

# take one observation from each strata
data.test <- data.table(Melanoma)[, .SD[1], by = c("invasion", "ici")]
setkeyv(data.test, c("invasion","ici"))

# identify the last event time for each strata
epsilon <- min(diff(unique(fit.coxph$y[,"time"])))/10
baseHazStrata <- as.data.table(predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE))
dt.times <- baseHazStrata[, .(beforeLastTime = time[.N]-epsilon, LastTime = time[.N], afterLastTime = time[.N]+epsilon), by = strata]

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
    
    prediction <- predict(fit.CSC, times = test.times, newdata = data.test, cause = 1)
    expect_equal(unname(is.na(prediction[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction[Ttempo,])), c(FALSE, FALSE, TRUE))
  }
})

test_that("Prediction with Cox model (strata) - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predictCox(fit.coxph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.double(prediction$hazard), 0)
  expect_equal(as.double(prediction$cumHazard), 0)
  expect_equal(as.double(prediction$survival), 1)
  
  prediction <- predictCox(fit.cph, times = test.times, newdata = Melanoma[1,,drop = FALSE])
  expect_equal(as.double(prediction$hazard), 0)
  expect_equal(as.double(prediction$cumHazard), 0)
  expect_equal(as.double(prediction$survival), 1)
})

test_that("Prediction with CSC (strata)  - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction), 0)
})

#### 3- [predictCSC] survtype = "survival" ####
cat("predictCSC, survtype = \"survival\" \n")
set.seed(10)
d <- sampleData(3e2, outcome = "competing.risks")
d$time <- round(d$time,2)
ttt <- sample(x = unique(sort(d$time)), size = 10) 
d2 <- d
times <- sort(c(0,d$time))

test_that("Prediction with CSC (survtype = survival)  - no strata",{
  CSC.fit <- CSC(Hist(time,event)~ X1+X2+X6,data=d, method = "breslow", survtype = "survival")
  
  p1 <- predict(CSC.fit, newdata = d2, times = times, cause = 1)
  pGS <- pec::predictEventProb(CSC.fit, newdata = d2, times = times, cause = 1)
  
  # expect_equal(unname(p1), unname(pGS))
})

test_that("Prediction with CSC (survtype = survival)  - strata",{
  CSC.fitS <- CSC(Hist(time,event)~ strata(X1) + X5 + strata(X3) + X7 +X2,data=d, method = "breslow", survtype = "survival")
  
  p1 <- predict(CSC.fitS, newdata = d2, times = times, cause = 1)
  pGS <- pec::predictEventProb(CSC.fitS, newdata = d2, times = times, cause = 1)
  
  # expect_equal(max(abs(na.omit(p1 - pGS)))>1e-8,FALSE)
})



#### 4- [predictCox] Dealing with weights ####
cat("weigthed Cox model \n")
set.seed(10)
data(Melanoma)
wdata <- runif(nrow(Melanoma), 0, 1)
times1 <- unique(Melanoma$time)

fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma)
fitW.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, weights = wdata)

fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fitW.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE, weights = wdata)

# res <- predictCox(fit.coxph, times = Melanoma$time, newdata = Melanoma)
test_that("Prediction with Cox model - weights",{
expect_error(resW <- predictCox(fitW.coxph, times = Melanoma$time, newdata = Melanoma))
expect_error(resW <- predictCox(fitW.cph, times = Melanoma$time, newdata = Melanoma))
})
# resGS <- survival:::predict.coxph(fit.coxph, times = times1, newdata = Melanoma, type = "expected")
# resGSW <- survival:::predict.coxph(fitW.coxph, times = times1, newdata = Melanoma, type = "expected")

# expect_equal(diag(res$cumHazard), resGS)
# expect_equal(diag(resW$cumHazard), resGSW)


#### 5- [predictCox,CSC] Check influence of the order of the prediction times ####
cat("Order of the prediction times \n")
data(Melanoma)
times2 <- sort(c(0,0.9*min(Melanoma$time),Melanoma$time[5],max(Melanoma$time)*1.1))
newOrder <- sample.int(length(times2),length(times2),replace = FALSE)

test_that("Prediction with Cox model - sorted vs. unsorted times",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick, data = Melanoma)
  predictionUNS <- predictCox(fit.coxph, times = times2[newOrder], newdata = Melanoma[1:5,], keep.times = FALSE)
  predictionS <- predictCox(fit.coxph, times = times2, newdata = Melanoma[1:5,], keep.times = FALSE)
  # predictSurvProb(fit.coxph, times = times2[newOrder], newdata = Melanoma)
  # predictSurvProb(fit.coxph, times = times2, newdata = Melanoma)
  expect_equal(predictionS, lapply(predictionUNS, function(x){x[,order(newOrder)]}))

  fit.cph <- cph(Surv(time,status == 1) ~ thick, data = Melanoma, y = TRUE, x = TRUE)
  predictionUNS <- predictCox(fit.cph, times = times2[newOrder], newdata = Melanoma[1:5,], keep.times = FALSE)
  predictionS <- predictCox(fit.cph, times = times2, newdata = Melanoma[1:5,], keep.times = FALSE)
  expect_equal(predictionS, lapply(predictionUNS, function(x){x[,order(newOrder)]}))
})

test_that("Prediction with CSC - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick, data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1, keep.times = FALSE)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1, keep.times = FALSE)
  expect_equal(predictionS, predictionUNS[,order(newOrder)])
})

test_that("Prediction with Cox model (strata) - sorted vs. unsorted times",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion), data = Melanoma)
  predictionUNS <- predictCox(fit.coxph, times = times2[newOrder], newdata = Melanoma, keep.times = FALSE, keep.strata = FALSE)
  predictionS <- predictCox(fit.coxph, times = times2, newdata = Melanoma, keep.times = FALSE, keep.strata = FALSE)
  expect_equal(predictionS, lapply(predictionUNS, function(x){x[,order(newOrder)]}))
  
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion), data = Melanoma, y = TRUE, x = TRUE)
  predictionUNS <- predictCox(fit.cph, times = times2[newOrder], newdata = Melanoma, keep.times = FALSE, keep.strata = FALSE)
  predictionS <- predictCox(fit.cph, times = times2, newdata = Melanoma, keep.times = FALSE, keep.strata = FALSE)
  expect_equal(predictionS$hazard, predictionUNS$hazard[,order(newOrder)])
  expect_equal(predictionS$cumHazard, predictionUNS$cumHazard[,order(newOrder)])
  expect_equal(predictionS$survival, predictionUNS$survival[,order(newOrder)])
})
# na.omit(predictionS$hazard)

test_that("Prediction with CSC (strata) - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion), data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1)
  expect_equal(predictionS, predictionUNS[,order(newOrder)])
})

test_that("Deal with negative time points",{
  expect_equal(unname(predictCox(fit.coxph, times = -1, newdata = dataset1)$survival), matrix(1,nrow = nrow(dataset1), ncol = 1))
  expect_equal(unname(predict(fit.CSC, times = -1, newdata = dataset1, cause = 1)), matrix(0,nrow = nrow(dataset1), ncol = 1))
})

test_that("Deal with NA in times",{
  expect_error(predictionS <- predictCox(fit.coxph, times = c(times2,NA), newdata = Melanoma))
  expect_error(predictionS <- predict(fit.CSC, times = c(times2,NA), newdata = Melanoma, cause = 1))
})


#### 6- [predictCox] Check baseline hazard ####
cat("Estimation of the baseline hazard \n")
data(Melanoma)

test_that("baseline hazard - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumHazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumHazard, 
               basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumHazard[predictCox(fit.cph)$hazard>0], 
               basehaz(fit.cph)$hazard, tolerance = 1e-8)
  
  ## possible differences due to different fit - coef(fit.coxph)-coef(fit.cph)
  expect_equal(predictCox(fit.cph),
               predictCox(fit.coxph, centered = TRUE), 
               tolerance = 100*max(abs(coef(fit.coxph)-coef(fit.cph))))
})

#### strata
test_that("baseline hazard (strata) - order of the results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, keep.strata = TRUE)$strata, 
               basehaz(fit.coxph)$strata)
  expect_equal(predictCox(fit.cph, keep.strata = TRUE)$strata[predictCox(fit.cph)$hazard>0], 
               basehaz(fit.cph)$strata)
  # expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), 
  #              as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata))
})

test_that("baseline hazard (strata) - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumHazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumHazard[predictCox(fit.cph)$hazard>0], 
               basehaz(fit.cph)$hazard)
  
  ## !!! not the same ordering in the strata between cph and coxph thus the results in predictCox differ in order
  ## possible differences due to different fit - coef(fit.coxph)-coef(fit.cph)
  # expect_equal(basehaz(fit.coxph)$hazard[predictCox(fit.cph)$hazard>0], 
  #              basehaz(fit.cph)$hazard)
  # expect_equal(predictCox(fit.cph, centered = TRUE),
  #              predictCox(fit.coxph, centered = TRUE), 
  #              tolerance = 100*max(abs(coef(fit.coxph)-coef(fit.cph))))
})

#### 7- [predictCox] Check format of the output ####
cat("Format of the output \n")
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

## strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)

test_that("baseline hazard (strata) - order of the results",{
  expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), as.numeric(basehaz(fit.coxph)$strata))
  expect_equal(as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata)[predictCox(fit.cph, keep.strata = TRUE)$hazard>0], as.numeric(basehaz(fit.cph)$strata))
})

test_that("baseline hazard (strata) - correct number of events",{ 
  strata <- interaction(Melanoma$invasion, Melanoma$ici)
  timePerStrata <- tapply(fit.coxph$y[,"time"],strata, function(x){length(unique(x))})
  
  lengthRes <- unlist(lapply(predictCox(fit.coxph, keep.times = TRUE, keep.strata = TRUE), length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
  lengthRes <- unlist(lapply(predictCox(fit.cph, keep.times = TRUE, keep.strata = TRUE), length))
  expect_equal(unname(lengthRes), rep(sum(timePerStrata), 5))
})

test_that("Prediction with Cox model (strata) - export of strata and times",{
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

test_that("Prediction with Cox model (strata) - consistency of hazard/cumHazard/survival",{
  predictTempo <- predictCox(fit.coxph, times = times1, newdata = dataset1)
  expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumHazard,1,diff)), tolerance = 1e-8)
  expect_equal(predictTempo$survival, exp(-predictTempo$cumHazard), tolerance = 1e-8)
})

predictTempo <- predictCox(fit.coxph, times = c(0,times1[1:10]), newdata = dataset1[1:2,])
expect_equal(predictTempo$hazard[,-1], t(apply(predictTempo$cumHazard,1,diff)), tolerance = 1e-8)
expect_equal(predictTempo$survival, exp(-predictTempo$cumHazard), tolerance = 1e-8)


test_that("Prediction with Cox model (strata) - incorrect strata",{
  dataset1$invasion <- "5616"
  expect_error(predictCox(fit.coxph, times = times1, newdata = dataset1))
})

#### 9- Conditional CIF ####
cat("Conditional CIF \n")

d <- SimCompRisk(1e2)
d$time <- round(d$time,1)
ttt <- sample(x = unique(sort(d$time)), size = 10)
d2 <- SimCompRisk(1e2)

#### coxph function
CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow" )

test_that("Conditional CIF identical to CIF before first event", {
  pred <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, t0 = min(d$time)-1e-1)
  expect_equal(pred, predC)

  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, t0 = min(d$time)-1e-1)
  expect_equal(pred, predC)
})

test_that("Conditional CIF is NA after the last event", {
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, t0 = max(d$time)+1)
  expect_equal(all(is.na(predC)), TRUE)
  
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, t0 = max(d$time)+1)
  expect_equal(all(is.na(predC)), TRUE)
  
  t0 <- mean(range(d$time))
  ttt0 <- c(t0,ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt0, t0 = t0)
  expect_equal(all(is.na(predC[,ttt0<t0])), TRUE)
  expect_equal(all(!is.na(predC[,ttt0>=t0])), TRUE)
})

test_that("Value of the conditional CIF", {
  ttt <- sort(c(0,ttt))
  indexT0 <- 5
  
  cumH1 <- predictCox(CSC.fit$models$`Cause 1`, newdata = d2, times = ttt[indexT0]-1e-6)[["cumHazard"]]
  cumH2 <- predictCox(CSC.fit$models$`Cause 2`, newdata = d2, times = ttt[indexT0]-1e-6)[["cumHazard"]]
  Sall <- exp(-cumH1-cumH2)
  
  predRef <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt[indexT0]-1e-6)
  
  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt)
  predC_manuel <- (pred-as.double(predRef))/as.double(Sall)
  predC_manuel[,seq(1,indexT0-1)] <- NA
  
  predC_auto <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, t0 = ttt[indexT0])
  expect_equal(predC_auto,predC_manuel)
})


#### 10- Others ####
cat("Others \n")
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

test_that("Prediction with CSC - categorical cause",{
predictRisk(CSC.h, newdata = dn[1,], times = c(2), cause = "1")
})

predictRisk(CSC.h, newdata = dn, times = c(1,2,3.24,3.25,3.26,5,10,15,20), cause = cause)
predictRisk(CSC.s, newdata = dn, times = c(1,5,10,15,20), cause = cause)

predictCox(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20))

predictRisk(CSC.h$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[1]], newdata = dn, times = c(5,10,15,20), cause = cause)

predictRisk(CSC.h$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)
predictRisk(CSC.s$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause)


#----------------------------------------------------------------------
### predictRisk.R ends here
