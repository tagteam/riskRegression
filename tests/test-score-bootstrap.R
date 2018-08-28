library(riskRegression)
library(survival)
library(testthat)
test_that("loob binary",{
    learndat=sampleData(200,outcome="binary")
    lr1a = glm(Y~X6,data=learndat,family=binomial)
    lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
    ## leave-one-out bootstrap
    set.seed(5)
    loob.se0 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=FALSE)
    set.seed(5)
    loob.se1 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=TRUE)
    expect_equal(loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
    expect_equal(loob.se0$Brier$contrasts$delta,loob.se1$Brier$contrasts$delta)
})

test_that("bootcv binary (multi.state.test)",{
    learndat=sampleData(200,outcome="binary")
    lr1a = glm(Y~X6,data=learndat,family=binomial)
    lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
    ## leave-one-out bootstrap
    set.seed(5)
    bootcv.se0 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se1 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se2 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=TRUE)
    set.seed(5)
    bootcv.se3 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=TRUE)
    bootcv <- list(bootcv.se0,bootcv.se1,bootcv.se2,bootcv.se3)
    ## delta
    for (i in 1:4)
        for (j in 2:4) 
            for (m in c("AUC","Brier"))
                expect_equal(bootcv[[i]][[m]]$contrasts$delta,bootcv[[j]][[m]]$contrasts$delta)
    ## lower, upper
    for (m in c("AUC","Brier"))
        expect_equal(bootcv[[2]][[m]]$contrasts[,.(lower,upper)],bootcv[[4]][[m]]$contrasts[,.(lower,upper)])
    ## p-value (does not work yet)
    ## for (m in c("AUC","Brier")){
        ## expect_equal(bootcv[[3]][[m]]$contrasts[,.(p)],bootcv[[4]][[m]]$contrasts[,.(p)])
    ## }
})

if(FALSE){ ## [:failed test:]
    test_that("loob survival",{
        learndat=sampleData(200,outcome="survival")
        cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
        cox2a = coxph(Surv(time,event)~X7+X8+X9,data=learndat,x=TRUE,y=TRUE)
        ## leave-one-out bootstrap
        set.seed(5)
        loob.se0 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=FALSE)
        set.seed(5)
        loob.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=TRUE)
        loob.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="loob",B=100,se.fit=TRUE,metrics="brier",conservative=TRUE)
        expect_equal(loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
        expect_equal(loob.se0$Brier$contrasts$delta,loob.se1$Brier$contrasts$delta)
    })
}

if(FALSE){ ## [:failed test:]
test_that("bootcv survival (multi.state.test)",{
    learndat=sampleData(200,outcome="survival")
    learndat[,eventtime:=NULL]
    learndat[,censtime:=NULL]
    cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
    cox2a = coxph(Surv(time,event)~X7+X8+X9,data=learndat,x=TRUE,y=TRUE)
    ## rf1 = rfsrc(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=learndat,ntree=500)
    ## rf1 = rfsrc(Surv(time,event)~.,data=learndat,ntree=500)
    ## leave-one-out bootstrap
    set.seed(5)
    bootcv.se0 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se1 <- Score(list("COX1"=cox1a,"COX2"=cox2a,"RF"=rf1),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=FALSE)
    set.seed(5)
    bootcv.se2 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=TRUE,conservative=TRUE)
    set.seed(5)
    bootcv.se3 <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=TRUE,conservative=TRUE)
    bootcv <- list(bootcv.se0,bootcv.se1,bootcv.se2,bootcv.se3)
    ## delta
    for (i in 1:4)
        for (j in 2:4) 
            for (m in c("AUC","Brier"))
                expect_equal(bootcv[[i]][[m]]$contrasts$delta,bootcv[[j]][[m]]$contrasts$delta)
    ## lower, upper
    for (m in c("AUC","Brier"))
        expect_equal(bootcv[[2]][[m]]$contrasts[,.(lower,upper)],bootcv[[4]][[m]]$contrasts[,.(lower,upper)])
    ## p-value
    ## for (m in c("AUC","Brier"))
        ## expect_equal(bootcv[[3]][[m]]$contrasts[,.(p)],bootcv[[4]][[m]]$contrasts[,.(p)])
})
}
