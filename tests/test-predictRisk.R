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

set.seed(10)
n <- 300
df <- SimCompRisk(n)
dn <- SimCompRisk(17)
seqTime <- c(unique(sort(df$time)), max(df$time) + 1)
for(model in c("coxph","cph")){
    for(method.ties in c("breslow","efron")){
        CSC.NS <- CSC(Hist(time,event) ~ X1*X2,data = df, method = method.ties, fitter = model)
        for(cause in 1:2){ 
            Event0.NS <- predictRisk(CSC.NS, newdata = dn, times = seqTime, cause = cause)
            EventTest.NS <- predict(CSC.NS, newdata = dn, times = seqTime, cause = cause)
            test_that(paste0("pec vs predictRiskRR - no strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
                expect_equal(as.numeric(Event0.NS),as.numeric(unlist(EventTest.NS)), tolerance = 1e-8)
            })
        }
    }
}

#### 2- With strata ####
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)
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
            test_that(paste0("pec vs predictRisk - strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
                expect_equal(as.numeric(Event0.S),as.numeric(unlist(EventTest.S)), tolerance = 1e-8)
            })
        }
    }
}


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


library(riskRegression)
library(data.table)
library(rbenchmark)
library(prodlim)
library(rms)
library(testthat)
library(survival)
library(pec)
#library(riskRegression)


#### 1- Without strata ####
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
        test_that(paste0("predictSurvProb vs predictRisk - without strata, method.ties/method.baseHaz/model = ",method.ties,"/",model,""),{
            expect_equal(object = Surv0, expected = 1-SurvProb, tolerance=1e-8)
        })
    }
}

#### 2- With strata ####
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
        test_that(paste0("predictSurvProb vs predictRisk - without strata, method.ties/model = ",method.ties,"/",model,""),{
            a <- Surv0
            b <- 1-SurvProb
            attributes(a) <- attributes(b)
            expect_equal(object = a, expected = b, tolerance=1e-6)
        })
    }
}
SurvProb <- pec::predictSurvProb(coxph, ttt, newdata = df[17:18,])

#----------------------------------------------------------------------
### predictRisk.R ends here
