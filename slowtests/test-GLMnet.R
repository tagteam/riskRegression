### test-GLMnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: okt 17 2025 (08:17) 
## Version: 
## Last-Updated: mar 11 2026 (12:52) 
##           By: Thomas Alexander Gerds
##     Update #: 8
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
library(riskRegression)
library(prodlim)
test_that("penalized logistic regression", {
    m = lvm(X2~X1+X8)
    distribution(m,~X2) = binomial.lvm()
    regression(m,X2~X1+X8) = c(-.05,.07)
    set.seed(9)
    d = setDT(sim(m,1100))
    fit = glm(X2~X1+X8,data = d,family = "binomial")
    lasso_0 = GLMnet(X2~X1+X8,data = d,family = "binomial",alpha = 1,lambda = 0)
    lasso = GLMnet(X2~X1+X8,data = d,family = "binomial",alpha = 1)
    ridge = GLMnet(X2~X1+X8,data = d,family = "binomial",alpha = 0,cv = TRUE)
    expect_output(print(lasso))
    expect_output(print(ridge$penalty.factor.used))
})

test_that("penalized logistic regression with rcs, unpenalized, and different penalty parameter values", {
    set.seed(9)
    d <- sampleData(174)
    fit = glm(X2~X1+X8,data = d,family = "binomial")
    un = GLMnet(X2~unpenalized(X1)+X8,data = d,family = "binomial",alpha = 1)
    expect_equal(un$penalty.factor,c(X8 = 1,X11 = 0))
    splun = GLMnet(X2~unpenalized(X1)+pen(X3,4)+rcs(X8)+rcs(X9,4),data = d,family = "binomial",alpha = 1)
    expect_equal(splun$penalty.factor,
                 c(X31 = 4, X11 = 0, rcs_X8_term1 = 2, rcs_X8_term2 = 2, rcs_X9_term1 = 2, 
                   rcs_X9_term2 = 2, rcs_X9_term3 = 2))
    ridge2 = GLMnet(X2~pen(X1,0)+pen(X8,2),data = d,family = "binomial",alpha = 0)
    expect_equal(ridge2$penalty.factor,c(X11 = 0, X8 = 2))
    predictRisk(splun,newdata = sampleData(8))
})


######################################################################
### test-GLMnet.R ends here
