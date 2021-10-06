library(testthat)
library(riskRegression)
library(timereg)
library(prodlim)
library(survival)
context("Binomial regression")
data(Melanoma, package = "riskRegression")
# {{{ "Absolute risk regression"
test_that("Absolute risk regression",{
    set.seed(17)
    d <- sampleData(300,outcome="competing.risks")
    f1 <- ARR(Hist(time,event)~X1+X2+strata(X6),data=d,cause=1)
    f2 <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2)+X6,data=d,cause=1,model="rcif")
    f3 <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2)+X6,data=d,cause=1,model="rcif2")
    f4 <- ARR(Hist(time,event)~X1+X2+strata(X8),data=d,cause=1)    
    f5 <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2)+X8,data=d,cause=1,model="rcif")
    f6 <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2)+X8,data=d,cause=1,model="rcif2")
    expect_equal(as.numeric(f1$timeConstantEffects$coef),c(f3$gamma))
    expect_equal(as.numeric(f4$timeConstantEffects$coef),c(f6$gamma))
    d[,X4:=factor(X4)]
    system.time(A <- ARR(Hist(time,event)~X1+X3+X4,data=d,cause=1))
    ## system.time(B <- timereg::comp.risk(Event(time,event)~ const(factor(X1))+ const(factor(X3)),data=d,cause=1,model="rcif"))
    system.time(B <- timereg::comp.risk(Event(time,event)~ const(X1)+ const(X3)+const(X4),data=d,cause=1,model="rcif2"))
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})

# }}}
# {{{ "Logistic risk regression"
test_that("Logistic risk regression",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- LRR(Hist(time,event)~X1+X2,data=d,cause=1)
    b <- timereg::comp.risk(Event(time,event)~ const(X1)+const(X2),data=d,cause=1,model="logistic2")
    expect_equal(as.numeric(a$timeConstantEffects$coef),c(b$gamma))
    A <- LRR(Hist(time,event)~X1+strata(X2),data=d,cause=1)
    B <- timereg::comp.risk(Event(time,event)~ const(X1)+X2,data=d,cause=1,model="logistic2")
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})

# }}}
# {{{ "Censoring model"
test_that("Censoring model",{
    ## NOTE: expect_error with regexp=NA tests no error 
    expect_error(ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cause=1), regexp = NA) 
    expect_error(ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~thick+strat(invasion)+epicel,cause=1), regexp = NA)
    expect_error(ARR(Hist(time,status)~sex+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~logthick+age,cause=1), regexp = NA)
})
# }}}
