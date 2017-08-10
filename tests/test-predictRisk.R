### test-predictRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 10 2017 (08:56) 
## Version: 
## Last-Updated: Aug 10 2017 (11:21) 
##           By: Thomas Alexander Gerds
##     Update #: 3
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

test_that("Additional arguments: example with imputation of missing data", {
    library(randomForestSRC)
    library(survival)
    data(pbc,package="survival")
    set.seed(10)
    forest <- rfsrc(Surv(time,status)~chol+age+sex,data=pbc,ntree=10,nsplit=10)
    ## no specification
    expect_error(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km"))
    ## correct specification
    expect_output(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km",predictRisk.args=list("rfsrc"=list(na.action="na.impute")))))
    ## wrong specification
    expect_error(print(Score(list(forest),formula=Hist(time,status)~1,data=pbc,conf.int=FALSE,metrics="brier",censModel="km",predictRisk.args=list("randomForestSRC"=list(na.action="na.impute")))))
})

######################################################################
### test-predictRisk.R ends here
