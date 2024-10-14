### test-wglm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 21 2021 (11:04) 
## Version: 
## Last-Updated: mar 15 2023 (13:14) 
##           By: Brice Ozenne
##     Update #: 25
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Packages
library(testthat)
library(data.table)
library(mets)
library(riskRegression)
library(survival)

## * Simulate data
set.seed(10)
n <- 250
tau <- 1:5
d <- sampleData(n, outcome = "competing.risks")
dFull <- d[event!=0] ## remove censoring
dSurv <- d[event!=2] ## remove competing risk

## * no censoring
test_that("wglm - no censoring",{
    test <- wglm(regressor.event = ~ X1 +  X8, formula.censor = Surv(time,event==0) ~ 1,
                 times = tau, data = dFull)

    GS <- suppressWarnings(logitIPCW(formula = Event(time,event)~X1 + X8,
                                     time = tau[5], data = dFull))

    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"], tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,summary(test, print = FALSE)[[5]][,"Std. Error"], summary(GS)$coef[,"Std.Err"], tolerance = 1e-5)
 
    ## ate
    test.ate <- ate(test, data = dFull, times = tau, treatment = "X1", verbose = FALSE)
    GS.ate <- logitATE(formula = Event(time,event)~X1 + X8,
                       time = tau[5], data = dFull, treat.model = X1~1)
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,estimate], unname(GS.ate$difriskG),
                 tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,estimate], unname(GS.ate$difriskG),
                 tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,se], unname(GS.ate$se.difriskG),
                 tolerance = 1e-5)
})

## * right censoring (but no competing risks)
test_that("wglm - censoring",{
    #### no covariate in censoring model ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ 1,
                 times = tau, data = dSurv, product.limit = TRUE)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event)~X1 + X8,
                                     cens.model = ~1,
                                     time = tau[5], data = dSurv, cens.code = 0, cause = 1))
    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"], tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,summary(test, print = FALSE)[[5]][,"Std. Error"], summary(GS)$coef[,"Std.Err"], tolerance = 1e-5)

    ## check weights
    e.KM <- as.data.table(predictCoxPL(coxph(Surv(time,event==0) ~ 1, data = dSurv),
                                       type = c("hazard","survival")))

    ## ate
    test.ate <- ate(test, data = dSurv, times = tau, treatment = "X1", verbose = FALSE)
    GS.ate <- suppressWarnings(logitATE(formula = Event(time,event)~X1 + X8, cens.model = ~1, cens.code = 0, cause = 1,
                                        time = tau[5], data = dSurv, treat.model = X1~1, binreg = FALSE))
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,estimate], unname(GS.ate$difriskG),
                 tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,se], unname(GS.ate$se.difriskG),
                 tolerance = 1e-5)

    #### stratified censoring model ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ X1,
                 times = tau, data = dSurv, product.limit = FALSE)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event) ~ X1 + X8,
                                     cens.model = ~X1,
                                     time = tau[5], data = dSurv, cens.code = 0, cause = 1))

    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"], tolerance = 1e-5)

    #### censoring model with continuous covariate ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ X1+X8,
                 times = tau, data = dSurv, product.limit = FALSE)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event) ~ X1 + X8,
                    cens.model = ~X1+X8,
                    time = tau[5], data = dSurv, cens.code = 0, cause = 1))

    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"], tolerance = 1e-5)
})

## * competing risks
test_that("wglm - competing risks",{
    #### no covariate in censoring model ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ 1,
                 times = tau, data = d, product.limit = TRUE, cause = 1)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event) ~ X1 + X8,
                    cens.model = ~1,
                    time = tau[5], data = d, cens.code = 0, cause = 1))
    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"], tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,summary(test, print = FALSE)[[5]][,"Std. Error"], summary(GS)$coef[,"Std.Err"], tolerance = 1e-5)

    ## ate
    test.ate <- ate(test, data = d, times = tau, treatment = "X1", verbose = FALSE)
    GS.ate <- suppressWarnings(logitATE(formula = Event(time,event)~X1 + X8, cens.model = ~1, cens.code = 0, cause = 1,
                                        time = tau[5], data = d, treat.model = X1~1, binreg = FALSE))
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,estimate], unname(GS.ate$difriskG),
                 tolerance = 1e-5)
    expect_equal(ignore_attr=TRUE,test.ate$diffRisk[5,se], unname(GS.ate$se.difriskG),
                 tolerance = 1e-5)

    #### stratified censoring model ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ X1,
                 times = tau, data = d, product.limit = FALSE, cause = 1)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event) ~ X1 + X8,
                    cens.model = ~X1,
                    time = tau[5], data = d, cens.code = 0, cause = 1))

    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"],
                 tolerance = 1e-5)
    ## expect_equal(summary(test, print = FALSE)[[5]][,"Std. Error"], summary(GS)$coef[,"Std.Err"], tolerance = 1e-5)

    #### censoring model with continuous covariate ####
    test <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ X1+X8,
                 times = tau, data = d, product.limit = FALSE)
    GS <- suppressWarnings(logitIPCW(formula = Event(time,event) ~ X1 + X8,
                    cens.model = ~X1+X8,
                    time = tau[5], data = d, cens.code = 0, cause = 1))

    expect_equal(ignore_attr=TRUE,coef(test, time = tau[5]), summary(GS)$coef[,"Estimate"],
                 tolerance = 1e-5)
    ## expect_equal(summary(test, print = FALSE)[[5]][,"Std. Error"], summary(GS)$coef[,"Std.Err"], tolerance = 1e-5)
})



## if(FALSE){ ## simulation study comparing standard errors - very slow!! (2.5 min with 10 cpus)
##     set.seed(10)
##     n <- 100
##     tau <- 3

##     warper <- function(i, n, tau){
##         iData <- sampleData(n, outcome = "competing.risks")

##         iModelRR <- wglm(regressor.event = ~ X1 + X8, formula.censor = Surv(time,event==0) ~ X1+X8,
##                          times = tau, data = iData, product.limit = FALSE)
##         iRR <- ate(iModelRR, treatment = "X1", data = iData, times = tau, cause = 1, verbose = FALSE)
##         iMets <- logitATE(formula = Event(time,event) ~ X1 + X8,
##                           cens.model = ~X1+X8, treat.model = X1~1,
##                           time = tau, data = iData, cens.code = 0, cause = 1)

##         iOut <- rbind(data.frame(method = "riskRegression", sim = i, n = n, estimate = confint(iRR)$diffRisk$estimate, se = confint(iRR)$diffRisk$se),
##                       data.frame(method = "mets", sim = i, n = n, estimate = iMets$difriskG, se = iMets$se.difriskG))
##         return(iOut)
##     }

##     library(pbapply)
##     ls.sim <- pblapply(1:1000, function(i){warper(i=1, n= n, tau = tau)}, cl = 10)
##     dt.sim <- as.data.table(do.call(rbind,ls.sim))

##     dt.sim[,.(rep = .N, se = sd(estimate), see = mean(se)),by="method"]
##     ##            method  rep        se       see
##     ## 1: riskRegression 1000 0.1092620 0.1036145
##     ## 2:           mets 1000 0.1092611 0.1039674
##     ## mets slight better with n=250 and tau=3
## }
##----------------------------------------------------------------------
### test-wglm.R ends here
