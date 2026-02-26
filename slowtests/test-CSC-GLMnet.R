### test-coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 28 2025 (09:31) 
## Version: 
## Last-Updated: feb 26 2026 (14:36) 
##           By: Thomas Alexander Gerds
##     Update #: 18
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
library(survival)
library(glmnet)
test_that("unpenalized Cox models",{
    set.seed(17)
    d <- sampleData(1000)
    a <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d)
    b <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0,lambda = 0)
    u <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0)
    expect_equal(as.numeric(coef(a$models[[1]])),as.numeric(b$models[[1]]$selected.beta),tolerance = 0.001)
    expect_true(all(abs(as.numeric(u$models[[1]]$selected.beta))<abs(as.numeric(b$models[[1]]$selected.beta))))
})

test_that("penalized Cox models",{
    set.seed(17)
    d <- sampleData(100)
    test <- sampleData(1000)
    # unpenalized
    a <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d)
    # ridge
    b <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0)
    # elnet
    e <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0.5)
    # lasso
    l <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 1)
    x = Score(list(unpenalized = a,rigde = b, elnet = e,lasso = l),
              formula = Hist(time,event)~1,
              data = test,
              summary = "risk",
              times = 3)
    riskRegression::predictRisk(l,newdata = test,times = 3)
    ## plotRisk(x,times = 3)
    ## plotRisk(x,times = 3,models = c("lasso","elnet"))
    suppressMessages(expect_output(print(x)))
})

test_that("Cox models with block penaly",{
    set.seed(17)
    d <- sampleData(100)
    test <- sampleData(1000)
    # unpenalized
    a <- CSC(list(Hist(time,event)~unpenalized(X1)+X6+pen(X7,3)+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = c("glmnet","coxph"))
    ## b <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = c("glmnet","coxph"),
             ## fitter_arguments = list(list(penalty.factor = ),NULL)
    x = Score(list(a = a,b = b),
              formula = Hist(time,event)~1,
              data = test,
              summary = "risk",
              times = 3)
    riskRegression::predictRisk(l,newdata = test,times = 3)
    ## plotRisk(x,times = 3)
    ## plotRisk(x,times = 3,models = c("lasso","elnet"))
    suppressMessages(expect_output(print(x)))
})


######################################################################
### test-coxnet.R ends here
