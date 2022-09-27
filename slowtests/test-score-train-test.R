#AUC with censoring covariates example (in the train-test situation)
library(riskRegression)
library(testthat)
library(cmprsk)
test_that("Score is working with conservative/not option for Binary data ",{
  set.seed(18)
  trainCR.binary <- sampleData(200,outcome="binary") ## censoring does not depend on covariates
  testCR.binary <- sampleData(500,outcome="binary") ##
  glm1 = glm(Y~X1+X2+X7+X9,data=trainCR.binary,family="binomial")
  glm2 = glm(Y~X1+X2+X7,data=trainCR.binary,family="binomial")
  x<-Score(list("m1"=glm1,"m2"=glm2),
           formula=Y~1,data=testCR.binary,se.fit=1L,conservative = TRUE)
  y<-Score(list("m1"=glm1,"m2"=glm2),
           formula=Y~1,data=testCR.binary,se.fit=1L)
  # y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
  #          formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",cens.model="Hal9001",conservative = TRUE)
  expect_output(print(x))
  expect_output(print(y))## you can have look that all the SEs are pretty close
})

test_that("Score is working with conservative/not option for survival data",{
  set.seed(18)
  trainCR.surv <- sampleData(200,outcome="survival") ## censoring does not depend on covariates
  testCR.surv <- sampleData(500,outcome="survival") ##
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.surv,x=TRUE)
  cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.surv,x=TRUE)
  x<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Surv(time,event)~1,data=testCR.surv,se.fit=1L,times=c(4),conservative = TRUE)
  y<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Surv(time,event)~1,data=testCR.surv,se.fit=1L,times=c(4))
  z<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~X2+X7,data=testCR.surv,se.fit=1L,times=c(4),conservative = TRUE)
  w<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Hist(time,event)~X2+X7,data=testCR.surv,se.fit=1L,times=c(4))
  # y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
  #          formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",cens.model="Hal9001",conservative = TRUE)
  expect_output(print(x))
  expect_output(print(y))
  expect_output(print(z))
  expect_output(print(w)) ## you can have look that all the SEs are pretty close
})

test_that("Score is working with conservative/not option for competing risk data",{
  set.seed(18)
  trainCR.comprisk <- sampleData(200,outcome="competing.risks") ## censoring does not depend on covariates
  testCR.comprisk <- sampleData(500,outcome="competing.risks") ##
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
  x<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4),conservative = TRUE)
  y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~1,data=testCR.comprisk,se.fit=1L,times=c(4))
  z<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),conservative = TRUE)
  w<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
           formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4))
  # y<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X1+X2)"=csc2),
  #          formula=Hist(time,event)~X2+X7,data=testCR.comprisk,se.fit=1L,times=c(4),metrics="AUC",cens.model="Hal9001",conservative = TRUE)
  expect_output(print(x))
  expect_output(print(y))
  expect_output(print(z))
  expect_output(print(w)) ## you can have look that all the SEs are pretty close
})

