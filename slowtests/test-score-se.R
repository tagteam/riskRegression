#AUC with censoring covariates example
set.seed(18)
devtools::load_all()
library(riskRegression)
library(testthat)
library(cmprsk)
trainCR <- sampleData(200,outcome="competing.risks")
testCR <- sampleData(500,outcome="competing.risks")
# Cause-specific Cox regression
csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
csc2 = CSC(Hist(time,event)~X1+X2,data=trainCR)
csc3 = CSC(Hist(time,event)~X1,data=trainCR)
test_that("AUC, with covariates in censoring, competing risk",{
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="AUC")
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~X2+X3,data=testCR,se.fit=1L,times=c(4),metrics="AUC")
  expect_equal(y$AUC$score$se,z$AUC$score$se,tolerance = 1e-2)
})

#AUC with competing risks
set.seed(18)
trainCR.comprisk <- sampleData(200,outcome="competing.risks")
testCR.comprisk <- sampleData(5000,outcome="competing.risks")
csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
csc2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
csc3 = CSC(Hist(time,event)~X1,data=trainCR.comprisk)
trainCR.survival <- sampleData(200,outcome="survival")
testCR.survival <- sampleData(5000,outcome="survival")
cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.survival,x=TRUE)
cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.survival,x=TRUE)
cox3 = coxph(Surv(time,event)~X1,data=trainCR.survival,x=TRUE)
test_that("Brier score SE against old implementation, no covariates in censoring, competing risk",{
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="Brier")
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="Brier",old.ic.method = TRUE)
  expect_equal(y$Brier$score$se,z$Brier$score$se,tolerance = 1e-4)
})

test_that("Brier score SE against old implementation, no covariates in censoring, survival",{
  y<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2,"cox(X1)"=cox3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="Brier")
  z<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2,"cox(X1)"=cox3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="Brier",old.ic.method = TRUE)
  expect_equal(y$Brier$score$se,z$Brier$score$se,tolerance = 1e-4)
})

test_that("AUC SE against old implementation, no covariates in censoring, competing risk",{
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="AUC")
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2,"CSC(X1)"=csc3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="AUC",old.ic.method = TRUE)
  expect_equal(y$AUC$score$se,z$AUC$score$se)
})

test_that("AUC SE against old implementation, no covariates in censoring, survival",{
  y<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2,"cox(X1)"=cox3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="AUC")
  z<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2,"cox(X1)"=cox3),
           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(4),metrics="AUC",old.ic.method = TRUE)
  expect_equal(y$AUC$score$se,z$AUC$score$se)
})


