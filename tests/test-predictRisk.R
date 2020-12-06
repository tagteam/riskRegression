### test-predictRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 10 2017 (08:56) 
## Version: 
## Last-Updated: Dec  5 2020 (13:20) 
##           By: Thomas Alexander Gerds
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(data.table)
library(lava)

## * SuperLearner
cat("[predictRisk.SuperPredictor] \n")
test_that("wrap the SuperLearner", {
    if (!requireNamespace("SuperLearner",quietly=TRUE)){
        message("Package SuperLearner not installed. Skip this test.")
    }else{
        library(SuperLearner)
        d <- sampleData(139,outcome="binary")
        sl = SuperLearner(Y = d$Y,
                          X = d[,c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"),with=0L],
                          family = binomial(),
                          SL.library = "SL.glm")
    }
})

## * rfsrc
cat("[predictRisk.rfsrc]  \n")
test_that("Additional arguments: example with imputation of missing data", {
    if (!requireNamespace("randomForestSRC",quietly=TRUE)){
        message("Package randomForestSRC not installed. Skip this test.")
    }else{
        data(pbc,package="survival")
        set.seed(10)
        forest <- rfsrc(Surv(time,status)~chol+age+sex,data=pbc,ntree=10,nsplit=10)
        ## no specification
        expect_error(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",cens.model="km"))
        ## correct specification
        expect_output(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",cens.model="km",predictRisk.args=list("rfsrc"=list(na.action="na.impute")))))
        ## wrong specification
        expect_error(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",cens.model="km",predictRisk.args=list("randomForestSRC"=list(na.action="na.impute")))))
    }
})

## * predictCox/predictCSC - to be reorganized (there is no clear test in this section)
cat("[predictRisk.CauseSpecificCox] \n")
test_that("Prediction with CSC - categorical cause",{
    set.seed(10)
    n <- 300
    df.S <- prodlim::SimCompRisk(n)
    df.S$time <- round(df.S$time,2)
    df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
    method.ties <- "efron"
    cause <- 1
    n <- 3
    set.seed(3)
    dn <- prodlim::SimCompRisk(n)
    dn$time <- round(dn$time,2)
    dn$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
    CSC.h3 <- CSC(Hist(time,event) ~ X1 + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph")
    CSC.h1 <- CSC(Hist(time,event) ~ strat(X1) + X3 + X2, data = df.S, ties = method.ties, fitter = "cph")
    CSC.h <- CSC(Hist(time,event) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph")
    CSC.s <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties, fitter = "coxph")
    predictRisk(CSC.h1, newdata = dn, times = c(5,10,15,20), cause = cause)
    predictRisk(CSC.h3, newdata = dn, times = c(5,10,15,20), cause = cause)
    CSC.h0 <- CSC(Hist(time,event) ~ X1 + X3 + X2, data = df.S, ties = method.ties, fitter = "cph")
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
    expect_success(predictRisk(CSC.s$models[[2]], newdata = dn, times = c(5,10,15,20), cause = cause))
})

## * [predictRisk.glm] vs. lava
cat("[predictRisk.glm] vs. lava \n")

## ** data
n <- 100
set.seed(10)
dt <- sampleData(n, outcome="binary")

fit <- glm(formula = Y ~ X1+X2, data=dt, family = "binomial")

## ** iid
test_that("[predictRisk.glm] compare to lava",{

    e.RR <- predictRisk(fit, newdata = dt, iid = TRUE, average.iid = FALSE)
    e.RR.iid <- attr(e.RR,"iid")
    attr(e.RR,"iid") <- NULL

    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        return(expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1")))
    })
    ## check point estimate
    expect_equal(e.RR, predict(fit, newdata = dt, type = "response"))
    expect_equal(unname(e.RR), unname(e.lava$coef))

    ## check variance
    expect_equal(unname(crossprod(e.RR.iid)), unname(e.lava$vcov))
})

## ** average.iid 
test_that("[predictRisk.glm] compare to lava (average.iid, no factor)",{

    ## no factor
    e.RR0 <- predictRisk(fit, newdata = dt, iid = TRUE, average.iid = FALSE)
    
    e.RR <- predictRisk(fit, newdata = dt, iid = TRUE, average.iid = TRUE)
    e.RR.average.iid <- attr(e.RR,"average.iid")
    attr(e.RR,"average.iid") <- NULL

    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        eXb <- expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1"))
        return(list(risk = eXb))},
        average=TRUE)

    ## check point estimate
    expect_equal(unname(mean(e.RR)), unname(e.lava$coef))

    ## check iid
    expect_equal(unname(rowMeans(attr(e.RR0,"iid"))), unname(e.RR.average.iid[,1]))

    ## check variance
    expect_equal(unname(sum((e.RR.average.iid + (e.RR-mean(e.RR))/NROW(e.RR))^2)), unname(e.lava$vcov)[1,1])
})

## ** average.iid with factor
test_that("[predictGLM] compare to lava (average.iid, factor)",{

    factor <- TRUE
    attr(factor,"factor") <- matrix(1:NROW(dt), ncol = 1)
    
    e.RR0 <- predictRisk(fit, newdata = dt, iid = TRUE, average.iid = FALSE)
    
    e.RR <- predictRisk(fit, newdata = dt, average.iid = factor)
    e.RR.average.iid <- attr(e.RR,"average.iid")
    attr(e.RR,"average.iid") <- NULL
    e.RR <- e.RR * (1:NROW(dt))
    
    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        eXb <- expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1"))
        return(list(risk = eXb*(1:NROW(data))))},
        average=TRUE)

    ## check point estimate
    expect_equal(unname(mean(e.RR)), unname(e.lava$coef))

    ## check iid
    expect_equal(unname(rowMeans(rowMultiply_cpp(attr(e.RR0,"iid"), scale = 1:NROW(dt)))), unname(e.RR.average.iid[[1]])[,1])

    ## check variance
    expect_equal(unname(sum((e.RR.average.iid[[1]] + (e.RR-mean(e.RR))/NROW(e.RR))^2)), unname(e.lava$vcov)[1,1])
})


######################################################################
### test-predictRisk.R ends here
