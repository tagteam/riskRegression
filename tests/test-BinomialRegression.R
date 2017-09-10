library(testthat)
library(riskRegression)
library(timereg)
library(prodlim)
library(survival)
context("Binomial regression")
data(Melanoma)

test_that("Absolute risk regression",{
    set.seed(17)
    d <- sampleData(300,outcome="competing.risks")
    a <- ARR(Hist(time,event)~X1+X2,data=d,cause=1)
    b <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2),data=d,cause=1,model="rcif")
    expect_equal(as.numeric(a$timeConstantEffects$coef),c(b$gamma))
    d[,X4:=factor(X4)]
    system.time(A <- ARR(Hist(time,event)~X1+X3+X4,data=d,cause=1))
    ## system.time(B <- timereg::comp.risk(Event(time,event)~ const(factor(X1))+ const(factor(X3)),data=d,cause=1,model="rcif"))
    system.time(B <- timereg::comp.risk(Event(time,event)~ const(X1)+ const(X3)+const(X4),data=d,cause=1,model="rcif"))
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})

test_that("Logistic risk regression",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- LRR(Hist(time,event)~X1+X2,data=d,cause=1)
    b <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2),data=d,cause=1,model="logistic")
    expect_equal(as.numeric(a$timeConstantEffects$coef),c(b$gamma))
    A <- LRR(Hist(time,event)~X1+strata(X2),data=d,cause=1)
    B <- timereg::comp.risk(Event(time,event)~ const(X1)+X2,data=d,cause=1,model="logistic")
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})

test_that("Censoring model",{
    f1 <- ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cause=1)
    f1a <- ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~thick+strat(invasion)+epicel,cause=1)
    f2 <- ARR(Hist(time,status)~sex+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~logthick+age,cause=1)
})

