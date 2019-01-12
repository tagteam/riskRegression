### test-Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (14:30) 
## Version: 
## last-updated: Dec 12 2018 (07:34) 
##           By: Thomas Alexander Gerds
##     Update #: 114
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
library(rms)
library(pec)
library(riskRegression)
library(pROC)
library(data.table)
library(Daim)
context("riskRegression")
# {{{ "R squared/IPA"
test_that("R squared/IPA", { 
    set.seed(112)
    d <- sampleData(43,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    r1 <- rsquared(f1,newdata=d)
    r2 <- IPA(f2,newdata=d)
    full <- Score(list(f1=f1,f2=f2),formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
    expect_equal(r1$IPA.gain[1],full$Brier$score[model=="f1",IPA])
    expect_equal(r2$IPA[2],full$Brier$score[model=="f2",IPA])
})

# }}}
# {{{ "vcov AUC"
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
    test_that("vcov AUC",{
        set.seed(112)
        d <- sampleData(43,outcome="binary")
        f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
        f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
        test <- Score(list(f1,f2),keep="vcov",formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
        expect_equal(dim(test$AUC$vcov),c(2,2))
        ## survival
        set.seed(112)
        d <- sampleData(112,outcome="survival")
        f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x=TRUE,y=TRUE)
        f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x=TRUE,y=TRUE)
        test <- Score(list(a=f1,f2),times=c(5,7),keep="vcov",formula=Surv(time,event)~1,data=d,conf.int=TRUE,metrics=c("brier","auc"))
        expect_equal(dim(test$AUC$vcov),c(4,4))
    })
}
# }}}
# {{{ "binary outcome: robustness against order of data set"
test_that("binary outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="binary")
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
# }}}
# {{{ "survival outcome: robustness against order of data set"
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
    test_that("survival outcome: robustness against order of data set",{
        set.seed(112)
        d <- sampleData(43,outcome="survival")
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
}
# }}}
# {{{ "competing risks outcome: check against pec"
test_that("competing risks outcome: check against pec",{
    set.seed(112)
    d <- sampleData(43,outcome="competing.risks")
    nd <- sampleData(43,outcome="competing.risks")
    library(pec)
    f <- FGR(Hist(time,event)~X1+X6,data=d,cause=1)
    a <- pec::pec(list(f),data=nd,times=c(2,5),formula=Hist(time,event)~1,cens.model="marginal",exact=FALSE)
    b <- Score(list(FGR=f),data=nd,formula=Hist(time,event)~1,cens.model="km",se.fit=FALSE,times=c(2,5),metrics="brier")
    expect_equal(a$AppErr$Reference[-1],b$Brier$score[model=="Null model",Brier])
    expect_equal(a$AppErr$FGR[-1],b$Brier$score[model=="FGR",Brier])
})

# }}}
# {{{ "competing risks outcome: robustness against order of data set"
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
test_that("competing risks outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="competing.risks")
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
})}
# }}}
# {{{ "survival outcome: Brier Score pec vs Score"
test_that("survival outcome: Brier Score pec vs Score",{
    set.seed(112)
    d <- sampleData(43,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
    p1 <- pec(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),exact=FALSE,start=NULL)
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=FALSE,metrics="brier")
    expect_equal(p1$AppErr$coxph,s1$Brier$score[model=="coxph",Brier])
    expect_equal(p1$AppErr$coxph.1,s1$Brier$score[model=="coxph.1",Brier])
    expect_equal(p1$AppErr$Reference,s1$Brier$score[model=="Null model",Brier])
})
# }}}
# {{{ "survival outcome: matrix input"
test_that("survival outcome: matrix input",{
    set.seed(112)
    dtrain <- sampleData(43,outcome="survival")
    dtest <- sampleData(4,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=dtrain, x = TRUE, y = TRUE)
    f2 <- predictRisk(f1,newdata=dtest,times=c(3,5,10))
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=dtest,times=c(3,5,10),conf.int=FALSE,null.model=0L,metrics="brier")
    expect_equal(s1$Brier$score[model=="coxph",Brier],s1$Brier$score[model=="matrix",Brier])
})

# }}}
# {{{ "survival outcome,Brier Score, external prediction"
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
    b <- Score(ModelsR,formula = Surv(time,status)~X1+X2,data=dat,times= c(5),se.fit=FALSE)
    cbind(b$Brier$score[,Brier],as.vector(unlist(a$AppErr)))
    expect_equal(b$Brier$score[,Brier],as.vector(unlist(a$AppErr)))
})
# }}}
# {{{ "survival outcome uncensored"
test_that("survival outcome uncensored",{
    library(mlbench)
    library(survival)
    library(randomForestSRC)
    data(BostonHousing)
    d=BostonHousing;d$status=1
    cx=coxph(Surv(medv,status)~.,d,x=TRUE)
    rfx=rfsrc(Surv(medv,status)~.,d,ntree=10)
    Score(list(Cox=cx,forest=rfx),data=d,formula=Hist(medv,status)~1,se.fit=TRUE)
}
# }}}
# {{{ "binary outcome: Brier"
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
    expect_equal(s1$Brier$score$Brier,s3$Brier$score$Brier)
    expect_equal(s2$Brier$score$Brier,s3$Brier$score$Brier)
    expect_equal(S1$Brier,S2$Brier)
    expect_equal(S1$Brier$score$Brier,S3$Brier$score$Brier)
    expect_equal(S2$Brier$score$Brier,S3$Brier$score$Brier)
})
# }}}
# {{{ "binary outcome: AUC"
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
# }}}
# {{{ "Leave one out bootstrap: Number of models and time points"
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
    test_that("Number of models and time points", {
        library(pec)
        data(GBSG2)
        setDT(GBSG2)
        ## fit1 <- coxph(Surv(time, cens)~horTh+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE)
        ## fit2 <- coxph(Surv(time, cens)~strata(horTh)+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE)
        fit1 <- cph(Surv(time, cens)~horTh+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        fit2 <- cph(Surv(time, cens)~strat(horTh)+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        GBSG2.test <- GBSG2
        setorder(GBSG2.test,time,-cens)
        ## predictCox(fit1,newdata=GBSG2.test,times=1000)
        r1 <- Score(list(a=fit2),data=GBSG2.test,times=1000,formula=Surv(time,cens)~1,plots="cali")
        set.seed(11)
        R1 <- Score(list(a=fit2),data=GBSG2.test,times=1000,B=50,split.method="loob",formula=Surv(time,cens)~1,plots="cali")
        setorder(GBSG2,time,cens)
        ## setorder(GBSG2.test,age)
        GBSG2 <- 7
        r2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(100,500,2000,1000),formula=Surv(time,cens)~1,plots="cali")
        set.seed(11)
        R2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(1000),B=50,split.method="loob",formula=Surv(time,cens)~1,plots="cali")
        ## r1$Calibration$plotframe
        ## r2$Calibration$plotframe[times==1000&model=="a"]
        ## r3 <- pec(list(a=fit2,b=fit1),data=GBSG2.test,exact=FALSE,times=c(1000),formula=Surv(time,cens)~1)
        expect_equal(r1$Brier$score[model=="a"],r2$Brier$score[model=="a" & times==1000])
        ## expect_equal(r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    })
    # {{{ "Bootstrap cross validation
    test_that("Number of models and time points", {
        library(pec)
        data(GBSG2)
        setDT(GBSG2)
        ## fit1 <- coxph(Surv(time, cens)~horTh+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE)
        ## fit2 <- coxph(Surv(time, cens)~strata(horTh)+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE)
        fit1 <- cph(Surv(time, cens)~horTh+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        fit2 <- cph(Surv(time, cens)~strat(horTh)+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        GBSG2.test <- GBSG2
        setorder(GBSG2.test,time,-cens)
        ## predictCox(fit1,newdata=GBSG2.test,times=1000)
        r1 <- Score(list(a=fit2),data=GBSG2.test,times=1000,formula=Surv(time,cens)~1,plots="cali")
        set.seed(11)
        R1 <- Score(list(a=fit2),data=GBSG2.test,times=1000,B=50,split.method="bootcv",formula=Surv(time,cens)~1,plots="cali")
        setorder(GBSG2,time,cens)
        ## setorder(GBSG2.test,age)
        GBSG2 <- 7
        r2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(100,500,2000,1000),formula=Surv(time,cens)~1,plots="cali")
        set.seed(11)
        R2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(1000),B=50,split.method="bootcv",formula=Surv(time,cens)~1,plots="cali")
        ## r1$Calibration$plotframe
        ## r2$Calibration$plotframe[times==1000&model=="a"]
        ## r3 <- pec(list(a=fit2,b=fit1),data=GBSG2.test,exact=FALSE,times=c(1000),formula=Surv(time,cens)~1)
        expect_equal(r1$Brier$score[model=="a"],r2$Brier$score[model=="a" & times==1000])
        ## expect_equal(r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    })
    # }}}
    # {{{ "LOOB: Number of models and time points"
    test_that("LOOB: Number of models and time points", {   
        library(testthat)
        library(survival)
        library(rms)
        library(pec)
        library(riskRegression)
        library(pec)
        library(pROC)
        library(data.table)
        library(Daim)
        data(GBSG2)
        setDT(GBSG2)
        setorder(GBSG2,time,-cens,age)
        fit1 <- cph(Surv(time, cens)~horTh+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        fit2 <- cph(Surv(time, cens)~strat(horTh)+age+menostat+tsize+pnodes+progrec+estrec, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        fit3 <- cph(Surv(time, cens)~horTh, data = GBSG2, x = TRUE,y=TRUE,surv=TRUE)
        setorder(GBSG2,time,-cens)
        set.seed(138)
        a <- Score(list(fit1=fit1,fit2=fit2,fit3=fit3),cens.model = "cox",data = GBSG2,times = c(2000,1000),metric = "brier",null.model = FALSE,contrast = FALSE,B=40,split.method="looboot",formula = Surv(time, cens)~horTh,se.fit=0L)
        setorder(GBSG2,time,cens)
        set.seed(138)
        b <- Score(list(fit2=fit2,fit1=fit1),cens.model = "cox",data = GBSG2,times = c(1000,2000),metric = "brier",null.model = FALSE,contrast = FALSE,B=40,split.method="looboot",formula = Surv(time, cens)~horTh,se.fit=0L) 
        expect_equal(a$Brier$score[model=="fit1"&times==2000],b$Brier$score[model=="fit1"&times==2000])
        expect_equal(a$Brier$score[model=="fit1"&times==1000],b$Brier$score[model=="fit1"&times==1000])
        expect_equal(a$Brier$score[model=="fit2"&times==2000],b$Brier$score[model=="fit2"&times==2000])
        expect_equal(a$Brier$score[model=="fit2"&times==1000],b$Brier$score[model=="fit2"&times==1000])
        B <- 50
        setorder(GBSG2,time,-cens)
        set.seed(138)
        A <- Score(list(fit1=fit1,fit2=fit2,fit3=fit3),cens.model = "cox",data = GBSG2,times = c(1000,365.25*4),metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
        print(A$Brier$score[model=="fit1"&times==365.25*4])
        setorder(GBSG2,time,cens)
        set.seed(138)
        A1 <- Score(list(fit1=fit1),cens.model = "cox",data = GBSG2,times = 365.25*4,metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
        print(A1$Brier$score[model=="fit1"&times==365.25*4])
        setorder(GBSG2,age)
        set.seed(138)
        A2 <- Score(list(fit2=fit2,fit1=fit1),cens.model = "cox",data = GBSG2,times = c(365.25*4,300,1000),metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
        print(A2$Brier$score[model=="fit1"&times==365.25*4])
        expect_equal(A$Brier$score[model=="fit1"&times==365.25*4],
                     A2$Brier$score[model=="fit1"&times==365.25*4])
        expect_equal(A$Brier$score[model=="fit1"&times==1000],
                     A2$Brier$score[model=="fit1"&times==1000])
        expect_equal(A$Brier$score[model=="fit2"&times==1000],
                     A2$Brier$score[model=="fit2"&times==1000])
        expect_equal(A1$Brier$score[model=="fit1"&times==365.25*4],
                     A2$Brier$score[model=="fit1"&times==365.25*4])
        expect_equal(A1$Brier$score[model=="fit1"&times==365.25*4],
                     A$Brier$score[model=="fit1"&times==365.25*4])
    })
}
# }}}
# }}}
#----------------------------------------------------------------------
### test-Score.R ends here
