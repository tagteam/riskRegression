### test-predictRisk-TD.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep  2 2017 (08:01) 
## Version: 
## Last-Updated: Dec  5 2020 (13:20) 
##           By: Thomas Alexander Gerds
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(testthat)
library(timereg)
library(survival)
test_that("start stop Cox: basehaz vs timereg",{
    data(veteran)
    ## remove ties
    set.seed(17)
    veteran$jtime <- veteran$time+rnorm(NROW(veteran),mean=0,sd=0.001)
    cox <- coxph(Surv(jtime, status) ~ celltype+karno + age + trt, data= veteran,x=1)
    ca <- cox.aalen(Surv(jtime, status) ~ prop(celltype)+prop(karno) + prop(age) + prop(trt), data= veteran)
    bh.cox <- basehaz(cox,centered=FALSE)
    ## timereg's cum does not include censored times where the hazard does not change
    expect_equal(bh.cox[!duplicated(bh.cox$hazard),]$time,ca$cum[-1,"time"])
    expect_equal(bh.cox[!duplicated(bh.cox$hazard),"hazard"],ca$cum[-1,"(Intercept)"])
    test.times <- c(0,4,300,8,veteran$jtime[c(17,10)])
    expect_equal(1-predictRisk(cox,newdata=veteran[c(47,8,17,10),],times=test.times),
                 predict(ca,se.fit=FALSE,newdata=veteran[c(47,8,17,10),],times=test.times)$S0)
    ## time-dependent
    vet2 <- survSplit(Surv(jtime, status) ~., veteran, cut=c(60, 120), episode ="timegroup")
    coxTD <- coxph(Surv(tstart, jtime, status) ~ celltype+karno + age + trt, data= vet2,x=1)
    caTD <- cox.aalen(Surv(tstart, jtime, status) ~ prop(celltype)+prop(karno) + prop(age) + prop(trt), data= vet2)
    bh.coxTD <- basehaz(coxTD,centered=FALSE)
    ## timereg's cum does not include censored times where the hazard does not change
    expect_equal(bh.coxTD[!duplicated(bh.coxTD$hazard),"hazard"],caTD$cum[-1,"(Intercept)"])
    ## with(basehaz(coxTD,centered=FALSE),plot(time,hazard,type="s"))
    ## lines(caTD$cum[,"time"],caTD$cum[,"(Intercept)"],col=2,type="s")
    expect_gt(unlist(riskRegression:::predictRisk.coxphTD(coxTD,newdata=vet2[c(17),],landmark=60,times=48)),
              unlist(riskRegression:::predictRisk.coxphTD(coxTD,newdata=vet2[c(17),],landmark=60,times=18)))
})

test_that("start stop CauseSpecificCox",{
    library(data.table)
    set.seed(9)
    d1 <- sampleData(201,outcome="competing.risk")
    d1[,start:=0]
    fit1 <- CSC(list(Hist(time=time,event=event,entry=start)~X1+X4+log(X7)+X9),data=d1)
    ## expect_equal(a <- riskRegression:::predictRisk.CSCTD(fit1,newdata=d1[c(17:32),],landmark=2:5,times=1,cause=1),
    ## b <- riskRegression:::predictRisk.CauseSpecificCox(fit1,newdata=d1[c(17:32),],landmark=2:5,times=1,cause=1))
    d2 <- sampleData(d1[,sum(event==0)],outcome="competing.risk")
    d2[,start:=d1[event==0,time]]
    d2[,time:=time+d1[event==0,time]]
    d3 <- sampleData(sum(d2$event==0),outcome="competing.risk")
    d3[,start:=d2[event==0,time]]
    d3[,time:=time+d2[event==0,time]]
    d <- rbindlist(list(d1,d2,d3))
    ## d[start>time,.(start,time,event)]
    fit <- CSC(list(Hist(time=time,event=event,entry=start)~X1+X4+log(X7)+X9),data=d)
    expect_success(riskRegression:::predictRisk.CSCTD(fit,newdata=d[c(17:32),],landmark=2:5,times=1))
})

######################################################################
### test-predictRisk-TD.R ends here
