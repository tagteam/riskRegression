### test-ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 15 2018 (11:42) 
## Version: 
## Last-Updated: maj  6 2019 (14:43) 
##           By: Brice Ozenne
##     Update #: 88
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(survival)
library(data.table)
library(testthat)
library(ipw)
context("Ate robust checks")

## * survival case

## ** no censoring - manual computation
n <- 5e1
tau <- 1.5
## simulate data
set.seed(10)
dtS <- sampleData(n,outcome="survival")
dtS$event <- 1
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$X1 <- as.numeric(dtS$X1)-1

## fit
test_that("check vs. manual calculations", {
    e.S <- coxph(Surv(time, event) ~ X1 + X2 + X3,data = dtS,x = TRUE)
    e.T <- glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")) ## dtS$X1
    e.ateRobust <- ateRobust(data = dtS, times = tau,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1, se = TRUE,
                             type = "survival")
    ## point estimate
    dtS1 <- data.table(X1 = 1, X2 = dtS$X2, X3 = dtS$X3)
    pred.logit <- riskRegression:::predictGLM(e.T, newdata = dtS1, average.iid = FALSE)
    iid.logit <- attr(pred.logit, "iid")
    attr(pred.logit, "iid") <- NULL
    pred.Surv <- predictCox(e.S, newdata = dtS1, times = tau, type = "survival", iid = TRUE)
    iid.Surv <- -pred.Surv$survival.iid[,1,]
    dtS[, pi := as.double(pred.logit)]
    dtS[, r := as.double(1-pred.Surv$survival)]
    expect_equal(dtS[,mean(r)],e.ateRobust$ate.value["risk.1","Gformula"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.IPCW"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.AIPCW"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.AIPCW_knownNuisance"])
  

    ## standard error
    iid.Gformula1 <- (dtS$r - mean(dtS$r))/NROW(dtS)
    iid.Gformula2 <- colMeans(iid.Surv)
    iid.Gformula <- iid.Gformula1 + iid.Gformula2
    expect_equal(sqrt(sum(iid.Gformula^2)),e.ateRobust$ate.se["risk.1","Gformula"])
    iid.IPTW1 <- (dtS[,X1*Y/pi] - dtS[,mean(X1*Y/pi)])/NROW(dtS)
    iid.IPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*dtS$Y/dtS$pi^2) )
    iid.IPTW <- iid.IPTW1 + iid.IPTW2
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.IPCW"])
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.AIPCW"])
    iid.AIPTW1 <- (dtS[,X1*Y/pi + r*(1-X1/pi)] - dtS[,mean(X1*Y/pi + r*(1-X1/pi))])/NROW(dtS)
    iid.AIPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*(dtS$Y-dtS$r)/dtS$pi^2) )
    iid.AIPTW3 <- colMeans( riskRegression::colMultiply_cpp(iid.Surv, scale = (1-dtS$X1/dtS$pi)) )
    iid.AIPTW <- iid.AIPTW1 + iid.AIPTW2 + iid.AIPTW3
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.AIPCW_knownNuisance"])
    ## (sqrt(sum(iid.AIPTW^2)) - sqrt(sum(iid.AIPTW1^2))) / sqrt(sum(iid.AIPTW^2))

    ## compare to ipw package
    e.ateRobust <- ateRobust(data = dtS, times = tau,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ X2+X3, se = TRUE,
                             type = "survival")
    RR.weight <- as.matrix(1/e.ateRobust$weight.treatment)
    
    ww <- ipwpoint(exposure=X1,family="binomial",link="logit",
                   denominator=~X2+X3,
                   data=dtS)$ipw.weights
    expect_equal(ww, RR.weight[1:NROW(RR.weight) + dtS$X1 * NROW(RR.weight)])
})

## ** censoring - agreement with ate
set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
e.cox <- coxph(Surv(time, event) ~ X1 + X2 + X3,
               data = dtS,
               x = TRUE)

test_that("Agreement ate-ateRobust (survival)",{
    e.ate <- ate(e.cox, treatment = "X1", times = 3, data = dtS, se = TRUE)

    e.ateRobust <- ateRobust(data = dtS, times = 3,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             se = TRUE, product.limit = FALSE,
                             type = "survival")

    e.ateRobustPL <- ateRobust(data = dtS, times = 3,
                               formula.event = Surv(time,event) ~ X1 + X2 + X3,
                               formula.censor = Surv(time,event==0) ~ X1,
                               formula.treatment = X1 ~ 1,
                               se = TRUE, product.limit = TRUE,
                               type = "survival")
    
    ## ateRobust with product.limit = FALSE agree with ate
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk],e.ate$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobust$ate.se[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk.se],e.ate$riskComparison[,diff.se]))

    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.19447, 0.54518, 0.35071), 
                  "IPTW.IPCW" = c(0.15098, 0.71242, 0.56145), 
                  "AIPTW.IPCW_knownNuisance" = c(0.1463, 0.75177, 0.60547), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.1463, 0.75177, 0.60547), 
                  "IPTW.AIPCW" = c(0.15097, 0.70596, 0.55499), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.14629, 0.7453, 0.59901),
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.14629, 0.7453, 0.59901))
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)

    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.04724, 0.13442, 0.12753), 
                  "IPTW.IPCW" = c(0.05683, 0.14767, 0.15823), 
                  "AIPTW.IPCW_knownNuisance" = c(0.05662, 0.15129, 0.15937), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.05571, 0.14053, 0.14949), 
                  "IPTW.AIPCW" = c(0.05681, 0.14873, 0.15921), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.05661, 0.15154, 0.15965), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0557, 0.14153, 0.15041)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
})


## * competing risk case

## ** no censoring - manual computation
n <- 1e2
tau <- 1.5

## simulate data
set.seed(10)
dtS <- sampleData(n,outcome="competing.risks")
dtS <- dtS[event!=0]
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$X1 <- as.numeric(dtS$X1)-1

## fit
e.S <- CSC(Hist(time, event) ~ X1 + X2 + X3,
           data = dtS)
e.T <- glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")) ## dtS$X1

e.ateRobust <- ateRobust(data = dtS, times = tau,
                         formula.event = Hist(time,event) ~ X1 + X2 + X3,
                         formula.censor = Hist(time,event==0) ~ X1,
                         formula.treatment = X1 ~ 1, se = TRUE, cause = 1,
                         type = "competing.risks")

## predict
dtS1 <- data.table(X1 = 1, X2 = dtS$X2, X3 = dtS$X3)
X <- as.double(model.matrix(formula(e.T), dtS1))

dprob.deta <- X * exp(-X*coef(e.T)) / (1+ exp(-X*coef(e.T)))^2
pred.logit <- 1 / (1 + exp(-X*coef(e.T)))
iid.logit <-  cbind(dprob.deta) %*% matrix(lava::iid(e.T), nrow = 1)

## GS <- attr(riskRegression:::predictGLM(e.T, newdata = dtS1, average.iid = FALSE), "iid")
## range(iid.logit - GS)

pred.risk <- predict(e.S, newdata = dtS1, times = tau, cause = 1, iid = TRUE)
iid.risk <- pred.risk$absRisk.iid[,1,]

dtS[, pi := pred.logit]
dtS[, r := as.double(pred.risk$absRisk)]

test_that("check point estimate vs. manual calculations", {
    expect_equal(dtS[,mean(r)],e.ateRobust$ate.value["risk.1","Gformula"])

    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.IPCW"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.AIPCW"])

    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.AIPCW_knownNuisance"])
})

test_that("check se vs. manual calculations", {

    iid.Gformula1 <- (dtS$r - mean(dtS$r))/NROW(dtS)
    iid.Gformula2 <- colMeans(iid.risk)
    iid.Gformula <- iid.Gformula1 + iid.Gformula2
    expect_equal(sqrt(sum(iid.Gformula^2)),e.ateRobust$ate.se["risk.1","Gformula"])

    iid.IPTW1 <- (dtS[,X1*Y/pi] - dtS[,mean(X1*Y/pi)])/NROW(dtS)
    iid.IPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*dtS$Y/dtS$pi^2) )
    iid.IPTW <- iid.IPTW1 + iid.IPTW2
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.IPCW"])
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.AIPCW"])

    iid.AIPTW1 <- (dtS[,X1*Y/pi + r*(1-X1/pi)] - dtS[,mean(X1*Y/pi + r*(1-X1/pi))])/NROW(dtS)
    iid.AIPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*(dtS$Y-dtS$r)/dtS$pi^2) )
    iid.AIPTW3 <- colMeans( riskRegression::colMultiply_cpp(iid.risk, scale = (1-dtS$X1/dtS$pi)) )
    iid.AIPTW <- iid.AIPTW1 + iid.AIPTW2 + iid.AIPTW3
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.AIPCW_knownNuisance"])
})


## ** agreement with ate
set.seed(10)
n <- 1e2
dtS <- sampleData(n,outcome="competing.risks")
## dtS[,min(time),by = event]

e.CSC <- CSC(Hist(time, event) ~ X1 + X2 + X3,
             data = dtS, surv.type = "hazard")

test_that("Agreement ate-ateRobust (competing.risks)",{
    ## NOT POSSIBLE: predictRisk does not recognise the argument product.limit
    ## e.ate <- ate(e.CSC, treatment = "X1", times = 3, data = dtS, cause = 1,
                 ## se = TRUE, product.limit = FALSE) 
    e.atePL <- ate(e.CSC, treatment = "X1", times = 3, data = dtS, cause = 1,
                   se = TRUE)

    ## data0 <- copy(dtS)
    ## data0[, X1 := factor("0", level = levels(dtS$X1))]
    ## data1 <- copy(dtS)
    ## data1[, X1 := factor("1", level = levels(dtS$X1))]
    ## mean(predict(e.CSC, cause = 1, newdata = data0, times = 3, product.limit = TRUE)$absRisk)
    
    e.ateRobust <- ateRobust(data = dtS, times = 3, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             product.limit = FALSE)
    e.ateRobustPL <- ateRobust(data = dtS, times = 3, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             product.limit = TRUE)

    expect_equal(as.double(e.ateRobustPL$ate.value[,"Gformula"]),
                 c(e.atePL$meanRisk[,meanRisk],e.atePL$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobustPL$ate.se[,"Gformula"]),
                 c(e.atePL$meanRisk[,meanRisk.se],e.atePL$riskComparison[,diff.se]))

    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28427, 0.30878, 0.02452), 
                  "IPTW.IPCW" = c(0.25855, 0.5098, 0.25125), 
                  "AIPTW.IPCW_knownNuisance" = c(0.25868, 0.50849, 0.24981), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.25868, 0.50849, 0.24981), 
                  "IPTW.AIPCW" = c(0.25835, 0.50462, 0.24627), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.25848, 0.50331, 0.24483), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.25848, 0.50331, 0.24483)
                  )
    
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5) 
    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.04577, 0.15642, 0.15524), 
                  "IPTW.IPCW" = c(0.04652, 0.16131, 0.16789), 
                  "AIPTW.IPCW_knownNuisance" = c(0.04636, 0.16792, 0.17302), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.04644, 0.15882, 0.16504), 
                  "IPTW.AIPCW" = c(0.04648, 0.1619, 0.16844), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.04632, 0.16799, 0.1731), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0464, 0.15946, 0.16567)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  "IPTW.IPCW" = c(0.25857, 0.5099, 0.25133), 
                  "AIPTW.IPCW_knownNuisance" = c(0.2587, 0.50859, 0.24988), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.2587, 0.50859, 0.24988), 
                  "IPTW.AIPCW" = c(0.25837, 0.50473, 0.24636), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.2585, 0.50341, 0.24491), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.2585, 0.50341, 0.24491)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.04576, 0.15641, 0.15524), 
                  "IPTW.IPCW" = c(0.04653, 0.16134, 0.16792), 
                  "AIPTW.IPCW_knownNuisance" = c(0.04636, 0.16804, 0.17314), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.04644, 0.15886, 0.16509), 
                  "IPTW.AIPCW" = c(0.04648, 0.16193, 0.16847), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.04632, 0.1681, 0.17322), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0464, 0.15951, 0.16572)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
   
})
######################################################################
### test-ateRobust.R ends here
