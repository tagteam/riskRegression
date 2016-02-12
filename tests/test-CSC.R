# Still need to check/update
# 1. When ties are present in own predictEventProb2
# - idea: the ties are calculated in CSC for efron ties - own code is breslow
# 2. Update own basehaz and the 2 predict functions to efron ties


library(pec)
library(CoxBoost)
library(survival)
library(riskRegression)
library(prodlim)

#
# Test data - basehaz:
#

d <- SimSurv(100)
d$X3 <- sample(1:12,100, replace=TRUE)
d$X4 <- sample(0:1,100, replace=TRUE)
time <- d$time
status <- d$status
Z <- d[,c("X1","X2","X3","X4")]

#
# Comparing function: 
#

f <- coxph(Surv(time,status)~X1+strata(X4),data=d, ties="breslow")

#
# Testing: own function against basehaz
#

all.equal(h_theta_breslow_strata(f,time,status,Z,cause=1)[,1:ncol(h_theta_breslow_strata(f,time,status,Z,cause=1))],
          basehaz(f, centered=FALSE)[,1:ncol(basehaz(f, centered=FALSE))])


#
# Testing predictEventProb
#


# Data: 
train <- SimCompRisk(100)
test <- SimCompRisk(10)
time <- train$time
status <- train$event
Z <- train[c("X1","X2")]
time2 <- round(train$time)

## cb.fit <- coxboost(Hist(time,cause)~X1+X2,cause=1,data=train,stepno=10)
## predictEventProb(cb.fit,newdata=test,times=seq(1:10),cause=1)

#with strate
cox.fit2  <- CSC(list(Hist(time,cause)~X2+strata(X1),Hist(time,cause)~X1+X2),data=train)
cox.fit3  <- CSC(Hist(time2,cause)~X1+strata(X2),data=train)

predictEventProb(cox.fit3,newdata=test,times=seq(0:10),cause=1)
predictEventProb2(cox.fit3,newdata=test,times=seq(0:10),cause=1,time,status,Z)

h_theta_breslow_strata(cox.fit3$models[[paste("Cause",1)]], time2, status, Z,cause=1)
basehaz(cox.fit3$models[[paste("Cause",1)]],centered=FALSE)

all.equal(unname(predictEventProb(cox.fit3,newdata=test,times=seq(0:10),cause=1)[,4]),
          PredictEventProb2(cox.fit3,newdata=test,times=seq(0:10),cause=1,time,status,Z)[,4])

