### test crossval/bootcv/loob for Score

library(devtools)
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
### Binary outcome 

## 5 fold cross-validation AUC and Brier with CI
test_that("test Score function with cross-validated AUC and Brier, B = 1",{
  set.seed(18)
  train = sampleData(400,outcome="binary")
  glm1 = glm(Y~X1+X2+X7+X9,data=train,family="binomial")
  glm2 = glm(Y~X3+X5+X6,data=train,family="binomial")
  x1 = Score(list("glm(X1+X2+X7+X9)"=glm1,"glm(X3+X5+X6)"=glm2),
             formula=Y~1,data=train,conf.int=TRUE,
             split.method="cv5",B=1)
  expect_output(print(x1))
})

test_that("test Score function with cross-validated AUC and Brier, B = 2",{
  set.seed(18)
  train = sampleData(400,outcome="binary")
  glm1 = glm(Y~X1+X2+X7+X9,data=train,family="binomial")
  glm2 = glm(Y~X3+X5+X6,data=train,family="binomial")
  x1 = Score(list("glm(X1+X2+X7+X9)"=glm1,"glm(X3+X5+X6)"=glm2),
             formula=Y~1,data=train,conf.int=TRUE,
             split.method="cv5",B=2)
  expect_output(print(x1))
})

# loob
test_that("test Score function with leave one out bootstrap AUC and Brier, B = 100",{
  set.seed(18)
  train = sampleData(4000,outcome="binary")
  glm1 = glm(Y~X1+X2+X7+X9,data=train,family="binomial")
  glm2 = glm(Y~X3+X5+X6,data=train,family="binomial")
  x1 = Score(list("glm(X1+X2+X7+X9)"=glm1,"glm(X3+X5+X6)"=glm2),
             formula=Y~1,data=train,conf.int=TRUE,
             split.method="loob",B=100)
  expect_output(print(x1))
})

# bootcv
test_that("test Score function with bootstrap crossvalidation AUC and Brier, B = 100",{
  set.seed(18)
  train = sampleData(400,outcome="binary")
  glm1 = glm(Y~X1+X2+X7+X9,data=train,family="binomial")
  glm2 = glm(Y~X3+X5+X6,data=train,family="binomial")
  x1 = Score(list("glm(X1+X2+X7+X9)"=glm1,"glm(X3+X5+X6)"=glm2),
             formula=Y~1,data=train,conf.int=TRUE,
             split.method="bootcv",B=100)
  expect_output(print(x1))
})

### Survival outcome

## 5 fold cross-validation AUC and Brier with CI
test_that("test Score function with cross-validated AUC and Brier, B = 1",{
  set.seed(18)
  trainSurv = sampleData(400,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
  cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
  x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
                            formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=5,
                            split.method="cv5",B=1)
  expect_output(print(x1))
})

test_that("test Score function with cross-validated AUC and Brier, B = 2",{
  set.seed(18)
  trainSurv = sampleData(400,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
  cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
  x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
             formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=5,
             split.method="cv5",B=2)
  expect_output(print(x1))
})

# loob
test_that("test Score function with leave one out bootstrap AUC and Brier, B = 100",{
  set.seed(18)
  trainSurv = sampleData(400,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
  cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
  x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
             formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=5,
             split.method="loob",B=100)
  expect_output(print(x1))
})

# bootcv
test_that("test Score function with bootstrap crossvalidation AUC and Brier, B = 100",{
  set.seed(18)
  trainSurv = sampleData(400,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
  cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
  x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
             formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=5,
             split.method="bootcv",B=100)
  expect_output(print(x1))
})

### Competing risk (warns about being under construction and Loglik not converging.)

## 5 fold cross-validation AUC and Brier with CI
test_that("test Score function with cross-validated AUC and Brier, B = 1",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  x1<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
            formula=Hist(time,event)~1,data=trainCR.comprisk,conf.int=TRUE,times=4,
            split.method="cv5",B=1)
  expect_output(print(x1))
})

test_that("test Score function with cross-validated AUC and Brier, B = 2",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  x1<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
            formula=Hist(time,event)~1,data=trainCR.comprisk,conf.int=TRUE,times=4,
            split.method="cv5",B=2)
  expect_output(print(x1))
})

# loob
test_that("test Score function with leave one out bootstrap AUC and Brier, B = 100",{
  set.seed(18)
  trainCR.comprisk <- sampleData(4000,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  x1<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
            formula=Hist(time,event)~1,data=trainCR.comprisk,conf.int=TRUE,times=4,
            split.method="loob",B=100)
  expect_output(print(x1))
})

# bootcv
test_that("test Score function with bootstrap crossvalidation AUC and Brier, B = 100",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  x1<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
            formula=Hist(time,event)~1,data=trainCR.comprisk,conf.int=TRUE,times=4,
            split.method="bootcv",B=100)
  expect_output(print(x1))
})

library(riskRegression)
library(testthat)

## problem with AUC and small sample sizes! (also problem for GLM!)
test_that("test that AUC is always >= 0.5 in small sample sizes", {
  set.seed(18)
  trainSurv = sampleData(35,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1,data=trainSurv, y=TRUE, x = TRUE)
  cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
  x1 = Score(list("Cox(X1)"=cox1,"Cox(X3+X5+X6)"=cox2),
             formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=seq(2,8,1),
             split.method="cv10",B=10,seed=9,metrics="AUC",summary="risk")
})
