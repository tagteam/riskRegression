### test-predictCox-CSC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:35) 
## Version: 
## last-updated: jun  3 2018 (23:09) 
##           By: Brice Ozenne
##     Update #: 151
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Setting
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(mstate)
tmat <- trans.comprisk(2, names = c("0", "1", "2"))

context("Risk prediction")
cat("prediction before and after the last event \n")

data(Melanoma)
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

#### no strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick*age, data = Melanoma, y = TRUE, x = TRUE)
fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma, fitter = "cph")


# {{{ "Prediction with CSC - NA after last event"
test_that("Prediction with CSC - NA after last event",{
  test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.vector(is.na(prediction$absRisk)), c(FALSE, FALSE, TRUE))
})

# }}}

# {{{ "Prediction with CSC - no event before prediction time"
test_that("Prediction with CSC - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction$absRisk), 0)
})

# }}}


#### strata
fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion) + strat(ici), data = Melanoma, fitter = "cph")

# take one observation from each strata
data.test <- data.table(Melanoma)[, .SD[1], by = c("invasion", "ici")]
setkeyv(data.test, c("invasion","ici"))

# identify the last event time for each strata
epsilon <- min(diff(unique(fit.coxph$y[,"time"])))/10
pred.coxph <- predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE)
baseHazStrata <- as.data.table(pred.coxph[c("times","hazard","cumhazard","strata","survival")])
dt.times <- baseHazStrata[, .(beforeLastTime = times[.N]-epsilon,
                              LastTime = times[.N],
                              afterLastTime = times[.N]+epsilon),
                          by = strata]

# }}}


# {{{ "Prediction with CSC (strata) - NA after last event"
test_that("Prediction with CSC (strata) - NA after last event",{
  for(Ttempo in 1:nrow(dt.times)){
    test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
    
    prediction <- predict(fit.CSC, times = test.times, newdata = data.test, cause = 1)
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
  }
})

# }}}


# {{{ "Prediction with CSC (strata)  - no event before prediction time"
test_that("Prediction with CSC (strata)  - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction$absRisk), 0)
})
# }}}


# {{{ 3- [predictCox,CSC] Check influence of the order of the prediction times

cat("Order of the prediction times \n")
data(Melanoma)
times2 <- sort(c(0,0.9*min(Melanoma$time),Melanoma$time[5],max(Melanoma$time)*1.1))
newOrder <- sample.int(length(times2),length(times2),replace = FALSE)


# {{{ "Prediction with CSC - sorted vs. unsorted times"
test_that("Prediction with CSC - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick, data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1, keep.times = FALSE)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1, keep.times = FALSE)
  expect_equal(predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
})

# }}}

# {{{ "Prediction with Cox model (strata) - sorted vs. unsorted times"
test_that("Prediction with Cox model (strata) - sorted vs. unsorted times",{
    times2 <- c(100,200,500)
    newOrder <- c(3,2,1)
    fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion), data = Melanoma, y = TRUE,  x = TRUE)
    predictionUNS <- predictCox(fit.coxph,times = c(500,200,100),newdata = Melanoma[1,],keep.times = FALSE,keep.strata = FALSE)
    predictionS <- predictCox(fit.coxph,times = c(100,200,500),newdata = Melanoma[1,],keep.times = FALSE,keep.strata = FALSE)
    expect_equal(predictionS[predictionS$type], lapply(predictionUNS[predictionUNS$type], function(x){x[,order(c(500,200,100)),drop=FALSE]}))
    fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion), data = Melanoma, y = TRUE, x = TRUE)
    predictionUNS <- predictCox(fit.cph,times = c(500,200,100),newdata = Melanoma[1,],keep.times = FALSE,keep.strata = FALSE)
    predictionS <- predictCox(fit.cph,times = c(100,200,500),newdata = Melanoma[1,],keep.times = FALSE,keep.strata = FALSE)
    expect_equal(predictionS[predictionS$type], lapply(predictionUNS[predictionUNS$type], function(x){x[,order(c(500,200,100)),drop=FALSE]}))
})
# na.omit(predictionS$hazard)

# }}}

# {{{ "Prediction with CSC (strata) - sorted vs. unsorted times"
test_that("Prediction with CSC (strata) - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion), data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1)
  expect_equal(predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
})

# }}}

# {{{ "Deal with negative time points"
test_that("Deal with negative time points",{
    expect_equal(unname(predictCox(fit.coxph, times = -1, newdata = dataset1)$survival),
                 matrix(1,nrow = nrow(dataset1), ncol = 1))
    expect_equal(unname(predict(fit.CSC, times = -1, newdata = dataset1, cause = 1)$absRisk),
                 matrix(0,nrow = nrow(dataset1), ncol = 1))
})
# }}}

# {{{ "Deal with NA in times"
test_that("Deal with NA in times",{
  expect_error(predictionS <- predictCox(fit.coxph, times = c(times2,NA), newdata = Melanoma))
  expect_error(predictionS <- predict(fit.CSC, times = c(times2,NA), newdata = Melanoma, cause = 1))
})

# }}}

# {{{ 4- Conditional CIF 
cat("Conditional CIF \n")

set.seed(10)
d <- SimCompRisk(1e2)
d$time <- round(d$time,1)
ttt <- sample(x = unique(sort(d$time)), size = 10)
d2 <- SimCompRisk(1e2)

#### coxph function
CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow")

# }}}

# {{{ "Conditional CIF identical to CIF before first event"
test_that("Conditional CIF identical to CIF before first event", {
  pred <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, landmark = min(d$time)-1e-5)
  expect_equal(pred, predC)

  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, landmark = min(d$time)-1e-5)
  expect_equal(pred, predC)
})

# }}}

# {{{ "Conditional CIF is NA after the last event"
test_that("Conditional CIF is NA after the last event", {
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, landmark = max(d$time)+1)
  expect_equal(all(is.na(predC$absRisk)), TRUE)
  
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, landmark = max(d$time)+1)
  expect_equal(all(is.na(predC$absRisk)), TRUE)
  
  t0 <- mean(range(d$time))
  ttt0 <- c(t0,ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt0, landmark = t0)
  expect_equal(all(is.na(predC$absRisk[,ttt0<t0])), TRUE)
  expect_equal(all(!is.na(predC$absRisk[,ttt0>=t0])), TRUE)
})

# }}}

# {{{ "Value of the conditional CIF | at the last event"
test_that("Value of the conditional CIF | at the last event", {
    tau <-  max(d[cause==2,time])
    
    predC_auto <- predict(CSC.fit, newdata = d2[1:5], cause = 2, times = tau, landmark = tau, product.limit = FALSE)
    pred <- as.data.table(predictCox(CSC.fit$models[[2]], times = tau, newdata = d2[1:5], type = "hazard"))
    expect_equal(as.double(predC_auto$absRisk),pred$hazard)

    predC_auto <- predict(CSC.fit, newdata = d2[1:5], cause = 2, times = tau, landmark = tau, product.limit = TRUE)
    pred <- as.data.table(predictCox(CSC.fit$models[[2]], times = tau, newdata = d2[1:5], type = "hazard"))
    expect_equal(as.double(predC_auto$absRisk),pred$hazard)
})

# }}}

# {{{ "Value of the conditional CIF"
test_that("Value of the conditional CIF", {
  sttt <- sort(c(0,ttt))
  indexT0 <- 5
    
  # product.limit = FALSE
  cumH1 <- predictCox(CSC.fit$models$`Cause 1`, newdata = d2, times = sttt[indexT0]-1e-6)[["cumhazard"]]
  cumH2 <- predictCox(CSC.fit$models$`Cause 2`, newdata = d2, times = sttt[indexT0]-1e-6)[["cumhazard"]]
  Sall <- exp(-cumH1-cumH2)
  
  predRef <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt[indexT0]-1e-6, product.limit = FALSE)
  
  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, product.limit = FALSE)
  predC_manuel <- (pred$absRisk-as.double(predRef$absRisk))/as.double(Sall)
  predC_manuel[,seq(1,indexT0-1)] <- NA
  
  predC_auto <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, landmark = sttt[indexT0], product.limit = FALSE)
  expect_equal(predC_auto$absRisk,predC_manuel)
  # predC_auto$absRisk - predC_manuel

  # product.limit = TRUE
  h1 <- predictCox(CSC.fit$models$`Cause 1`, newdata = d2, times = CSC.fit$eventTimes, type = "hazard")[["hazard"]]
  h2 <- predictCox(CSC.fit$models$`Cause 2`, newdata = d2, times = CSC.fit$eventTimes, type = "hazard")[["hazard"]]
  Sall <- apply(1-h1-h2,1, function(x){
      c(1,cumprod(x))[sindex(jump.times = CSC.fit$eventTimes, eval.times = sttt[indexT0]-1e-6)+1]
  })
  predRef <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt[indexT0]-1e-6, product.limit = TRUE)
  
  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, product.limit = TRUE)
  predC_manuel <- sweep(pred$absRisk-as.double(predRef$absRisk), MARGIN = 1, FUN ="/", STATS = as.double(Sall))
  predC_manuel[,seq(1,indexT0-1)] <- NA
  
  predC_auto <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, landmark = sttt[indexT0], product.limit = TRUE)
  expect_equal(predC_auto$absRisk,predC_manuel)
  # predC_auto$absRisk-predC_manuel
})
# }}}

# {{{ 6- Average influence function
cat("Average influence function \n")

## predictCox
set.seed(10)
d <- sampleData(80,outcome="competing.risks")
nd <- sampleData(4,outcome="competing.risks")

# }}}

# {{{ "average.iid - predict.CSC"
test_that("average.iid - predict.CSC", {

    ls.ff <- list(Hist(time,event)~X1 + X2 + X6,
                  Hist(time,event)~X1 + strata(X2) + X6)

    for(iF in ls.ff){ # iF <- ls.ff[[1]] 
        fit <- CSC(iF, data=d)
        
        resGS <- predict(fit, newdata=nd, times = 2:7, iid = TRUE, cause = 1)
        res1 <- predict(fit, newdata=nd, times = 2:7, average.iid = TRUE, store.iid = "full", cause = 1)
        res2 <- predict(fit, newdata=nd, times = 2:7, average.iid = TRUE, store.iid = "minimal", cause = 1)

        expect_equal(res1$absRisk.average.iid,res2$absRisk.average.iid)
        expect_equal(res2$absRisk.average.iid,t(apply(resGS$absRisk.iid,2:3,mean)))
    }

})

# }}}

## predictCSC
 
# }}}
# {{{ 7- absolute risk over 1
# this section does not perform any tests
# but show an example where the estimated absolute risk
# is over 1, probably because of the small sample size.
# I don't know if this is an issue.
if(FALSE){
    set.seed(5)
    d <- sampleData(80,outcome="comp")
    d[event==0,event:=1] # remove censoring
    ttt <- sort(unique(d[X1==1,time]))

    CSC.fit.s <- CSC(list(Hist(time,event)~ strata(X1)+X2+X9,
                          Hist(time,event)~ X2+strata(X4)+X8+X7),data=d)
    predict(CSC.fit.s,newdata=d[X1==1][1],cause=1,times=ttt,se=1L) # NA values due to absolute risk over 1
    predict(CSC.fit.s,newdata=d[X1==1][1],product.limit = FALSE,cause=1,times=ttt,se=1L) # NA values due to absolute risk over 1

    #### investigate how come we get absRisk > 1
    # since absRisk = int hazard1 * survival
    # it is possible if survival > 1
    #                or hazard(1) > 1
    
    ## restrict to the first strata and simplify the model
    ## to see if we can get an hazard > 1
    dt <- d[X1==1]
    setkeyv(dt, "time")
    m.coxph <- coxph(Surv(time,event>=1)~X2, data = dt, x = TRUE, ties = "breslow")

    # compute baseline hazard
    as.data.table(predictCox(m.coxph, centered = TRUE, type = "hazard")) # automatic
    rev(1/cumsum(rev(eXB))) # manual

    # normally they are normalized such that they are at most one:
    rev(rev(eXB)/cumsum(rev(eXB)))

    # but this does not work for new observations
    rev(unique(eXB)[1]/cumsum(rev(eXB))) # ok
    rev(unique(eXB)[2]/cumsum(rev(eXB))) # no hazard over 1

    # this also creates a problem when computing the suvival using the product limit estimator 
}
# }}}


#----------------------------------------------------------------------
### test-predictCox-CSC.R ends here
