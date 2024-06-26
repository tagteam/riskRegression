library(riskRegression)
library(survival)
library(data.table)
library(testthat)
library(randomForestSRC)

test_that("loob survival",{
    set.seed(8)
    ## learndat=sampleData(200,outcome="survival")
    learndat=sampleData(38,outcome="survival")
    cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
    cox2a = coxph(Surv(time,event)~X7+X8+X9,data=learndat,x=TRUE,y=TRUE)
    ## leave-one-out bootstrap
    set.seed(5)
    loob.se0 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=FALSE,progress.bar=NULL,verbose = -1)
    set.seed(5)
    loob.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=TRUE,progress.bar=NULL,verbose = -1)
    # covariate dependent censoring
    loob.se1a <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~X1,data=learndat,times=5,split.method="loob",B=100,se.fit=TRUE,progress.bar=NULL,verbose = -1)
    ## set.seed(5)
    ## loob.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=TRUE,metrics="brier",conservative=TRUE)
    expect_equal(ignore_attr=TRUE,loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
    expect_equal(ignore_attr=TRUE,loob.se0$Brier$contrasts$delta,loob.se1$Brier$contrasts$delta)
})

test_that("loob AUC comp risk",{
   
})

## GIVES WARNING
## Cannot do multi-split test with AUC yet. Forced multi.split.test=FALSE
test_that("bootcv survival (multi.state.test)",{
    set.seed(8)
    learndat=sampleData(200,outcome="survival")
    learndat[,eventtime:=NULL]
    learndat[,censtime:=NULL]
    cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
    cox2a = coxph(Surv(time,event)~X7+X8+X9,data=learndat,x=TRUE,y=TRUE)
    rf1 = rfsrc(Surv(time,event)~.,data=learndat,ntree=500)
    ## leave-one-out bootstrap
    set.seed(5)
    bootcv.se0 <- Score(list("COX1"=cox1a,"COX2"=cox2a,"RF"=rf1),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a,"RF"=rf1),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se2 <- Score(list("COX1"=cox1a,"COX2"=cox2a,"RF"=rf1),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=TRUE,conservative=TRUE)
    set.seed(5)
    bootcv.se3 <- Score(list("COX1"=cox1a,"COX2"=cox2a,"RF"=rf1),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=TRUE,conservative=TRUE)
    bootcv <- list(bootcv.se0,bootcv.se1,bootcv.se2,bootcv.se3)
    ## delta
    for (i in 1:4)
        for (j in 2:4)
            for (m in c("AUC","Brier"))
                expect_equal(ignore_attr=TRUE,bootcv[[i]][[m]]$contrasts$delta,bootcv[[j]][[m]]$contrasts$delta)
    ## lower, upper
    for (m in c("AUC","Brier"))
        expect_equal(ignore_attr=TRUE,bootcv[[2]][[m]]$contrasts[,.(lower,upper)],bootcv[[4]][[m]]$contrasts[,.(lower,upper)])
    ## p-value
    ## for (m in c("AUC","Brier"))
        ## expect_equal(ignore_attr=TRUE,bootcv[[3]][[m]]$contrasts[,.(p)],bootcv[[4]][[m]]$contrasts[,.(p)])
})

