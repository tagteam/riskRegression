library(riskRegression)
library(hal9001)
library(survival)
# warning these tests are very slow with continuous covariates (but one can also test with the continuous ones)
# all tests give a warning (as they should for now)

test_that("Predict risk is working with hal",{
    set.seed(18)
    trainCR.surv <- sampleData(270,outcome="survival") ## censoring does not depend on covariates
    testCR.surv <- sampleData(5,outcome="survival") ##
    cox1 = coxph(Surv(time,event)~X1+X9,data=trainCR.surv,y = TRUE,x = TRUE)
    hal1 = Hal9001(Surv(time,event)~X1+X9,data=trainCR.surv)
    p1 = predictRisk(hal1,newdata = testCR.surv[1,],times = c(1,3))
    p2 = predictRisk(cox1,newdata = testCR.surv,times = c(1,3))
    expect_output(print(p1))
})

# hal / train-test
test_that("Score is working with hal for survival data",{
    set.seed(18)
    trainCR.surv <- sampleData(200,outcome="survival") ## censoring does not depend on covariates
    testCR.surv <- sampleData(500,outcome="survival") ##
    cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.surv,x=TRUE)
    cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.surv,x=TRUE)
    suppressWarnings(x<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
                              formula=Surv(time,event)~X1+X2,
                              data=testCR.surv,
                              se.fit=1L,
                              times=c(4,5),
                              cens.model="Hal9001",
                              conservative = TRUE))
    expect_output(print(x))
})

test_that("Score is working with hal for competing risk data",{
  set.seed(18)
  trainCR.comprisk <- sampleData(200,outcome="competing.risks") ## censoring does not depend on covariates
  testCR.comprisk <- sampleData(500,outcome="competing.risks") ## censoring does not depend on covariates
  CSC1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  CSC2 = CSC(Hist(time,event)~X1+X2+X7,data=trainCR.comprisk)
  suppressWarnings(x<-Score(list("CSC(X1+X2+X7+X9)"=CSC1,"CSC(X1+X2)"=CSC2),
           formula=Hist(time,event)~X1+X2,data=testCR.comprisk,se.fit=1L,times=c(4),cens.model="Hal9001",conservative = TRUE))
  expect_output(print(x))
})

test_that("Score is working with hal for survival data / crossval",{ ## should also work with cvk and bootcv?
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="survival") ## censoring does not depend on covariates
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainCR.comprisk,x=TRUE)
  cox2 = coxph(Surv(time,event)~X1+X2+X7,data=trainCR.comprisk,x=TRUE)
  suppressWarnings(x<-Score(list("cox(X1+X2+X7+X9)"=cox1,"cox(X1+X2)"=cox2),
           formula=Surv(time,event)~X1+X2,data=trainCR.comprisk,conf.int=TRUE,times=4,
           split.method="loob",B=100,cens.model = "Hal9001",conservative = TRUE) )
  expect_output(print(x))
})


test_that("Score is working with hal for competing risk data / crossval",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  x<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
            formula=Hist(time,event)~X1+X2,data=trainCR.comprisk,conf.int=TRUE,times=4,
            split.method="loob",B=100,cens.model = "Hal9001",conservative = TRUE) 
  expect_output(print(x))
})

test_that("Score is working with hal for competing risk data / cv5, B = 1",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  suppressWarnings(x<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
           formula=Hist(time,event)~X1+X2,data=trainCR.comprisk,conf.int=TRUE,times=4,
           split.method="cv5",B=1,cens.model = "Hal9001",conservative = TRUE))
  expect_output(print(x))
})

test_that("Score is working with hal for competing risk data / cv5, B = 100",{
  set.seed(18)
  trainCR.comprisk <- sampleData(400,outcome="competing.risks")
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR.comprisk)
  csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR.comprisk)
  suppressWarnings(x<-Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2),
           formula=Hist(time,event)~X1+X2,data=trainCR.comprisk,conf.int=TRUE,times=4,
           split.method="cv5",B=1,cens.model = "Hal9001",conservative = TRUE) )
  expect_output(print(x))
})

