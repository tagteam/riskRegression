### test-predictSurv.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 28 2018 (10:06) 
## Version: 
## Last-Updated: sep 20 2019 (13:27) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(rms)
library(survival)
library(testthat)
context("predictSurv checks")

## * function predictSurv
## ** Simulate data
set.seed(10)
d <- sampleData(2e2, outcome = "competing.risk")
d.pred <- d[5:15]
seqTime <- c(0,d[["time"]][5:15],2.45,1e5)

## ** hazard case
e.CSC <- CSC(Hist(time, event)~ strata(X1) + strata(X2), data = d, surv.type = "hazard")
e.cox <- coxph(Surv(time, event>0)~ strata(X1) + strata(X2), data = d, x = TRUE, y = TRUE)
## e.cox <- coxph(Surv(time, event>0)~ 1, data = d, x = TRUE, y = TRUE)

test_that("predictSurv (hazard)", {
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = FALSE, iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime, iid = TRUE)
    expect_equal(GS$survival, test$survival)
    expect_equal(GS$survival.iid, test$survival.iid)

    testPL <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = TRUE, iid = TRUE)
    GSPL <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime, iid = TRUE)
    expect_equal(GSPL$survival, testPL$survival)
    expect_equal(GSPL$survival.iid, testPL$survival.iid)
})

## ** survival case
e.CSC <- CSC(Hist(time, event)~ X1 + X2, data = d, surv.type = "survival")
e.cox <- coxph(Surv(time, event>0)~ X1 + X2, data = d, x = TRUE, y = TRUE)

test_that("predictSurv (survival)", {
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = FALSE, iid = TRUE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime, iid = TRUE)
    expect_equal(GS$survival, test$survival)
    expect_equal(GS$survival.iid, test$survival.iid)

    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime, product.limit = TRUE, iid = TRUE)
    GS <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime, iid = TRUE)
    expect_equal(GS$survival, test$survival)
    expect_equal(GS$survival.iid, test$survival.iid)
})

## * function calcSurvival_cpp
prepareArgs <- function(object, data, times){

    jumpTimes <- object$eventTimes[object$eventTimes <= max(times)]
    
    out <- list()    
    out$ls_cumhazard <- lapply(object$models, function(iM){
        iH <- riskRegression::predictCox(iM, times = jumpTimes, keep.strata = TRUE, centered = FALSE)
        if(!is.null(iH$strata)){
            return(do.call(cbind,tapply(iH$cumhazard,iH$strata,list)))
        }else{
            return(cbind(iH$cumhazard))
        }
    })    
    out$eXb <- exp(cbind(riskRegression::coxLP(object$models[[1]], data = d.pred, center = FALSE),
                         riskRegression::coxLP(object$models[[2]], data = d.pred, center = FALSE)))
    out$nCause <- length(object$cause)
    out$theCause <- which(object$cause == object$theCause)-1
    out$hazardType <- (object$surv.type=="hazard")
    out$nNewObs <- NROW(data)
    out$nJumpTime <- length(jumpTimes)

    new.strata <- do.call(cbind,lapply(object$models, function(iM){ ## object$models <- 
        ls.infoVar <- riskRegression::coxVariableName(iM, model.frame = riskRegression::coxModelFrame(iM))
    
        riskRegression::coxStrata(iM, data = data,
                  sterms = ls.infoVar$strata.sterms, 
                  strata.vars = ls.infoVar$stratavars, 
                  strata.levels = ls.infoVar$strata.levels)
    }))
    out$Ustrata <- unique(new.strata)-1
    out$nStrata <- NROW(out$Ustrata)
    new.level.Ustrata <- apply(out$Ustrata+1,1,paste0,collapse="")
    new.Ustrata <- apply(new.strata,1,paste0,collapse="")
    out$ls_indexStrata <- lapply(new.level.Ustrata, function(iStrata){
        which(new.Ustrata==iStrata) - 1
    })
    return(out)
}

## ** Simulate data
set.seed(10)
d <- sampleData(75, outcome = "competing.risk")
d.pred <- d[5:15]
seqTime <- d[["time"]][5:15]
## seqTime <- c(0,d[["time"]][5:15],2.45,1e5)

## ** hazard case
expect_equal("[predictCSC] vs. calcSurvival_cpp - surv.type=\"hazard\"",{
    ## no strata
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "hazard")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime)
    test <- do.call(calcSurvBeforeJump_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime-(1e-12),
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)

    ## with strata
    e.CSC <- CSC(Hist(time, event)~ X6 + strata(X1), data = d, surv.type = "hazard")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime)
    test <- do.call(calcSurvBeforeJump_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime-(1e-12),
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)

    
    e.CSC <- CSC(list(Hist(time, event)~ X6 + strata(X1),
                      Hist(time, event)~ X6 + strata(X2)),
                 data = d, surv.type = "hazard")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime-(1e-12))
    test <- do.call(calcSurvBeforeJump_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime,
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)
})


## ** survival case
expect_equal("[predictCSC] vs. calcSurvival_cpp - surv.type=\"survival\"",{
    ## no strata
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime)
    test <- do.call(calcSurvival_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime,
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)

    ## with strata
    e.CSC <- CSC(Hist(time, event)~ X6 + strata(X1), data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime)
    test <- do.call(calcSurvival_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime,
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)

    
    e.CSC <- CSC(list(Hist(time, event)~ X6 + strata(X1),
                      Hist(time, event)~ X6 + strata(X2)),
                 data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    
    ls.args <- prepareArgs(e.CSC, data = d.pred, times = seqTime)
    test <- do.call(calcSurvival_cpp, ls.args)
    GS <- predict(e.CSC, type = "survival", times = jumpTime,
                  newdata = d.pred, product.limit = FALSE)
    expect_equal(GS$survival,test)
})


######################################################################
### test-predictSurv.R ends here
