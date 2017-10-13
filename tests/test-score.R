### test-Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (14:30) 
## Version: 
## last-updated: Oct 13 2017 (12:59) 
##           By: Thomas Alexander Gerds
##     Update #: 41
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
library(pec)
library(pROC)
library(data.table)
library(Daim)
context("riskRegression")

test_that("R squared", { 
    set.seed(112)
    d <- sampleData(112,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    full <- Score(list(f1,f2),formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
})

test_that("binary outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(112,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    f3 <- d$X8
    s1 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=TRUE)
    s1b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    setkey(d,X4)
    f3 <- d$X8
    s2 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95)
    s2b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    setorder(d,Y)
    f3 <- d$X8
    s3 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95)
    s3b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    ## lapply(names(s1),function(n){print(n);expect_equal(s1[[n]],s3[[n]])})
    s1$call$conf.int <- .95
    expect_equal(s1,s2)
    expect_equal(s1,s3)
    expect_equal(s1$AUC,s1b$AUC)
    expect_equal(s2$AUC,s2b$AUC)
    expect_equal(s3$AUC,s3b$AUC)
})
test_that("survival outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(112,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s1 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=TRUE)
    s1b <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    setkey(d,X4)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s2 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95)
    s2b <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    setorder(d,time,-event)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s3 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95)
    s3b <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    s1$call$conf.int <- .95
    expect_equal(s1,s2)
    expect_equal(s1,s3)
    expect_equal(s1$AUC,s1b$AUC)
    expect_equal(s2$AUC,s2b$AUC)
    expect_equal(s3$AUC,s3b$AUC)
})
test_that("competing risks outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(112,outcome="competing.risks")
    f1 <- CSC(Hist(time,event)~X1+X5+X8,data=d)
    f2 <- FGR(Hist(time,event)~X2+X6+X9+X10,data=d,cause=1)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s1 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=TRUE,cause=1)
    s1b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,cause=1,metrics="auc")
    setkey(d,X4)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s2 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,cause=1)
    s2b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,cause=1,metrics="auc")
    setorder(d,time,-event)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s3 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,cause=1)
    s3b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,cause=1,metrics="auc")
    s1$call$conf.int <- .95
    expect_equal(s1,s2)
    expect_equal(s1,s3)
    expect_equal(s1$AUC,s1b$AUC)
    expect_equal(s2$AUC,s2b$AUC)
    expect_equal(s3$AUC,s3b$AUC)
})
test_that("survival outcome: Brier Score pec vs Score",{
    set.seed(112)
    d <- sampleData(112,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
    p1 <- pec(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),exact=FALSE,start=NULL)
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=FALSE,metrics="brier")
    expect_equal(p1$AppErr$coxph,s1$Brier$score[model=="coxph",Brier])
    expect_equal(p1$AppErr$coxph.1,s1$Brier$score[model=="coxph.1",Brier])
    expect_equal(p1$AppErr$Reference,s1$Brier$score[model=="Null model",Brier])
})
test_that("survival outcome: matrix input",{
    set.seed(112)
    dtrain <- sampleData(112,outcome="survival")
    dtest <- sampleData(4,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=dtrain, x = TRUE, y = TRUE)
    f2 <- predictRisk(f1,newdata=dtest,times=c(3,5,10))
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=dtest,times=c(3,5,10),conf.int=FALSE,null.model=0L,metrics="brier")
    expect_equal(s1$Brier$score[model=="coxph",Brier],s1$Brier$score[model=="matrix",Brier])
})

test_that("survival outcome,Brier Score, external prediction",{
    ## generate simulated data
    set.seed(130971)
    n <- 4
    dat <- SimSurv(n)
    dat <- dat[order(dat$time,-dat$status),]
    ## define models
    ## Models <- list("Cox.X1" = coxph(Surv(time,status)~X1,data=dat,y=TRUE),
    Models <- list(
        constant=matrix(rep(0.43,n),ncol=1),
        "runif" = matrix(runif(n),ncol=1),"another.runif" = matrix(runif(n),ncol=1))
    ModelsR <- lapply(Models,function(x)1-x)
    ## training error
    a <- pec(Models,formula = Surv(time,status)~X1+X2,data=dat,times= c(5),exact=FALSE,start=NULL,verbose=TRUE)
    ## compare models
    b <- Score(ModelsR,formula = Surv(time,status)~X1+X2,data=dat,times= c(5))
    cbind(b$Brier$score[,Brier],as.vector(unlist(a$AppErr)))
    expect_equal(b$Brier$score[,Brier],as.vector(unlist(a$AppErr)))
})
test_that("binary outcome: Brier",{
    set.seed(47)
    D <- sampleData(n=47,outcome="binary")
    s1 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics="brier",cause="1")
    s2 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics=c("auc","brier"),cause="1")
    s3 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics=c("auc","brier"),se.fit=FALSE,cause="1")
    setkey(D,Y)
    S1 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics="brier",cause="1")
    S2 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics=c("auc","brier"),cause="1")
    S3 <- Score(list(X6=glm(Y~X6,data=D,family='binomial'),X9=glm(Y~X9,data=D,family='binomial'),X10=glm(Y~X10,data=D,family='binomial')),formula=Y~1,data=D,null.model=FALSE,metrics=c("auc","brier"),se.fit=FALSE,cause="1")
    expect_equal(s1,S1)
    expect_equal(s2,S2)
    expect_equal(s3,S3)
    expect_equal(s1$Brier,s2$Brier)
    expect_equal(s1$Brier,s3$Brier)
    expect_equal(s2$Brier,s3$Brier)
    expect_equal(S1$Brier,S2$Brier)
    expect_equal(S1$Brier,S3$Brier)
    expect_equal(S2$Brier,S3$Brier)
})

test_that("binary outcome: AUC", {   
    set.seed(17)
    y <- rbinom(100, 1, .5)
    x1 <- rnorm(100) + 1.5 * y
    x2 <- rnorm(100) + .5 * y
    x3 <- rnorm(100) + 2.5 * y
    x <- data.frame(x1,x2,x3)
    y <- as.factor(y)
    daimres <- Daim::deLong.test(x, labels=y, labpos="1")
    r1 <- pROC::roc(y~x1)
    r2 <- pROC::roc(y~x2)
    r3 <- pROC::roc(y~x3)
    procres <- pROC::roc.test(r1,r2)
    d <- data.frame(x1,x2,x3,y)
    ## Source(riskRegression)
    scoreres <- Score(list(X1=~x1,X2=~x2,X3=~x3),formula=y~1,data=d,null.model=FALSE,cause="1")
    ## Roc(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d)
    scoreres <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,cause="1")
    ## to avoid side effects of data.table features we check the following 
    scoreres1 <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,metrics="auc",cause="1")
    scoreres1a <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,metrics="auc",se.fit=0L,cause="1")
    expect_equal(scoreres$AUC,scoreres1$AUC)
    daim.auc <- daimres$AUC[,c("AUC","SD(DeLong)")]
    score.auc <- as.data.frame(scoreres$AUC$score[,c("AUC","se"),with=FALSE])
    rownames(score.auc) <- rownames(daim.auc)
    colnames(score.auc) <- colnames(daim.auc)
    expect_equal(daim.auc,score.auc)
    expect_equal(scoreres$AUC$score[["AUC"]],c(r1$auc,r2$auc,r3$auc))
    score.diff <- scoreres$AUC$contrasts[,c("delta.AUC","se","lower","upper","p"),with=FALSE]
    daim.diff <- daimres$difference
    expect_equal(daim.diff$"AUC Difference",-score.diff$delta.AUC)
    expect_equal(daim.diff$"CI(lower)",-score.diff$upper)
    expect_equal(daim.diff$"CI(upper)",-score.diff$lower)
    expect_equal(daim.diff$"P.Value",score.diff$p)
})

#----------------------------------------------------------------------
### test-Score.R ends here
