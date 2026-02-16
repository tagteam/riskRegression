### test-GLMnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: okt 17 2025 (08:17) 
## Version: 
## Last-Updated: okt 17 2025 (11:33) 
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
library(testthat)
library(lava)
library(data.table)
test_that("penalized logistic regression", {
    m = lvm(X2~X1+X8)
    distribution(m,~X2) = binomial.lvm()
    regression(m,X2~X1+X8) = c(-.05,.07)
    set.seed(9)
    d = setDT(sim(m,1100))
    d[,table(X2)]
    fit = glm(X2~X1+X8,data = d,family = "binomial")
    lasso = GLMnet(X2~X1+X8,data = d,family = "binomial",alpha = 1)
    ridge = GLMnet(X2~X1+X8,data = d,family = "binomial",alpha = 0)
    

})


######################################################################
### test-GLMnet.R ends here
