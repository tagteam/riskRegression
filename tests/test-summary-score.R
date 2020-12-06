### test-summary-score.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 12 2020 (07:48) 
## Version: 
## Last-Updated: Apr 12 2020 (09:46) 
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
library(survival)
library(riskRegression)
library(data.table)
context("riskRegression")
# {{{ "binary"
d1 <- sampleData(n=112,outcome="binary")
d2 <- sampleData(n=80,outcome="binary")
f1 <- glm(Y~X2+X8,data=d1,family="binomial")
f2 <- glm(Y~X1+X2+X5+X8+X6,data=d1,family="binomial")
# with null and se
x <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2)
xa <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="auc")
xb <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="brier")
print(x)
summary(x)
print(xa)
summary(xa)
print(xb)
summary(xb)
# without null with se
x <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,null.model=0L)
xa <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="auc",null.model=0L)
xb <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="brier",null.model=0L)
print(x)
summary(x)
print(xa)
summary(xa)
print(xb)
summary(xb)
# without null without se
x <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,null.model=0L,se.fit=0L)
xa <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="auc",null.model=0L,se.fit=0L)
xb <- Score(list(model.1=f1,model.2=f2),formula=Y~1,data=d2,metric="brier",null.model=0L,se.fit=0L)
print(x)
summary(x)
print(xa)
summary(xa)
print(xb)
summary(xb)



######################################################################
### test-summary-score.R ends here
