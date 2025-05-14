### test-coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 28 2025 (09:31) 
## Version: 
## Last-Updated: May 14 2025 (08:31) 
##           By: Thomas Alexander Gerds
##     Update #: 7
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
test_that("penalized Cox models",{
    set.seed(17)
    d <- sampleData(100)
    test <- sampleData(1000)
    a <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d)
    b <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0)
    e <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 0.5)
    l <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6+X9),data=d,fitter = "glmnet",
             alpha = 1)
    x = Score(list(unpenalized = a,rigde = b, elnet = e,lasso = l),
              formula = Hist(time,event)~1,
              data = test,
              summary = "risk",
              times = 3)
    ## plotRisk(x,times = 3)
    ## plotRisk(x,times = 3,models = c("lasso","elnet"))
    x
})


######################################################################
### test-coxnet.R ends here
