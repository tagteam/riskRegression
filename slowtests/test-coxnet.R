### test-coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 28 2025 (09:31) 
## Version: 
## Last-Updated: Apr 28 2025 (10:02) 
##           By: Thomas Alexander Gerds
##     Update #: 2
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
    a <- CSC(list(Hist(time,event)~X1+X6+X7+X8,Hist(time,event)~X1+X6),data=d)
    b <- CSC(Hist(time,event)~X1+X6+X7+X8+X9,data=d,fitter = "penalized")
    x = Score(list(unpenalized = a,penalize = b),formula = Hist(time,event)~1,data = d,summary = "risk",times = 3)
    plotRisk(x,times = 3)
})


######################################################################
### test-coxnet.R ends here
