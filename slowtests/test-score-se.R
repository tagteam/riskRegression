#AUC with censoring covariates example
devtools::load_all()
library(riskRegression)
library(testthat)
library(cmprsk)
test_that("AUC, covariates in censoring, competing risk, also an example of hal9001 implementation",{
  set.seed(18)
  trainCR.comprisk <- sampleData(200,outcome="competing.risks")
  testCR.comprisk <- sampleData(500,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
  x<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",conservative = TRUE)
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",cens.model="Hal9001",conservative = TRUE)
  expect_output(print(x))
})

#AUC with competing risks
# csc3 = CSC(Hist(time,event)~X1,data=trainCR.comprisk)
# cox3 = coxph(Surv(time,event)~X1,data=trainCR.survival,x=TRUE)
test_that("Brier score SE against old implementation, no covariates in censoring, competing risk",{
  set.seed(18)
  trainCR.comprisk <- sampleData(200,outcome="competing.risks")
  testCR.comprisk <- sampleData(5000,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="Brier")
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="Brier",old.ic.method = TRUE)
  expect_equal(y$Brier$score$Brier,z$Brier$score$Brier,tolerance = 1e-8)
  expect_equal(y$Brier$score$se,z$Brier$score$se,tolerance = 1e-4)
})

test_that("Brier score SE against old implementation, no covariates in censoring, survival",{
  set.seed(18)
  trainCR.survival <- sampleData(200,outcome="survival")
  testCR.survival <- sampleData(10000,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.survival,x=TRUE)
  cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.survival,x=TRUE)
system.time(y<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~1,data=testCR.survival,se.fit=1L,times=c(4),metrics="Brier"))
system.time( z<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~1,data=testCR.survival,se.fit=1L,times=c(4),metrics="Brier",old.ic.method = TRUE))
  expect_equal(y$Brier$score$se,z$Brier$score$se,tolerance = 1e-4)
})

test_that("AUC SE against old implementation, no covariates in censoring, competing risk",{
  set.seed(18)
  trainCR.comprisk <- sampleData(200,outcome="competing.risks")
  testCR.comprisk <- sampleData(5000,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC")
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",old.ic.method = TRUE)
  expect_equal(y$AUC$score$se,z$AUC$score$se,tolerance = 1e-5)
})

test_that("AUC SE against old implementation, no covariates in censoring, survival",{
  set.seed(18)
  trainCR.survival <- sampleData(200,outcome="survival")
  testCR.survival <- sampleData(5000,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.survival,x=TRUE)
  cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.survival,x=TRUE)
  y<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~1,data=testCR.survival,se.fit=1L,times=c(4),metrics="AUC")
  z<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~1,data=testCR.survival,se.fit=1L,times=c(4),metrics="AUC",old.ic.method = TRUE)
  expect_equal(y$AUC$score$se,z$AUC$score$se,tolerance = 1e-5)
})
