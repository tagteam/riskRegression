### test-coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 28 2025 (09:31) 
## Version: 
## Last-Updated: May 14 2025 (17:10) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
    x
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
    predictRisk(l,newdata = test,times = 3)
    ## plotRisk(x,times = 3)
    ## plotRisk(x,times = 3,models = c("lasso","elnet"))
    x
})


######################################################################
### test-coxnet.R ends here
