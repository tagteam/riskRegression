### test-summary-score.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 12 2020 (07:48) 
## Version: 
## Last-Updated: Jul 19 2021 (08:40) 
##           By: Thomas Alexander Gerds
##     Update #: 5
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
for (y in c("binary","survival","competing.risks")){
    print(y)
    d1 <- sampleData(n=112,outcome=y)
    d2 <- sampleData(n=80,outcome=y)
    if (y=="binary"){
        f1 <- glm(Y~X2+X8,data=d1,family="binomial")
        f2 <- glm(Y~X1+X2+X5+X8+X6,data=d1,family="binomial")
        ff <- Y~1
    }
    if (y=="survival"){
        f1 <- coxph(Surv(time,event)~X2+X8,data=d1,x=1L,y=1L)
        f2 <- coxph(Surv(time,event)~X1+X2+X5+X8+X6,data=d1,x=1L,y=1L)
        ff <- Hist(time,event)~1
    }
    if (y=="competing.risks"){
        f1 <- CSC(Hist(time,event)~X2+X8,data=d1)
        f2 <- CSC(Hist(time,event)~X1+X2+X5+X8+X6,data=d1)
        ff <- Hist(time,event)~1
    }
    # with null and se
    x <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2)
    xa <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="auc")
    xb <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="brier")
    print(x)
    summary(x)
    summary(x,what="contrast")
    summary(x,what="score")
    print(xa)
    summary(xa)
    summary(xa,what="contrast")
    summary(xa,what="score")
    print(xb)
    summary(xb)
    summary(xb,what="contrast")
    summary(xb,what="score")
    # without null with se
    x <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,null.model=0L)
    xa <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="auc",null.model=0L)
    xb <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="brier",null.model=0L)
    print(x)
    summary(x)
    summary(x,what="contrast")
    summary(x,what="score")
    print(xa)
    summary(xa)
    summary(xa,what="contrast")
    summary(xa,what="score")
    print(xb)
    summary(xb)
    summary(xb,what="contrast")
    summary(xb,what="score")
    # without null without se
    x <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,null.model=0L,se.fit=0L)
    xa <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="auc",null.model=0L,se.fit=0L)
    xb <- Score(list(model.1=f1,model.2=f2),formula=ff,data=d2,metric="brier",null.model=0L,se.fit=0L)
    print(x)
    summary(x)
    summary(x,what="contrast")
    summary(x,what="score")
    print(xa)
    summary(xa)
    summary(xa,what="contrast")
    summary(xa,what="score")
    print(xb)
    summary(xb)
    summary(xb,what="contrast")
    summary(xb,what="score")
}
# }}}


######################################################################
### test-summary-score.R ends here
