### test-Score_survival_outcome.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 17 2024 (12:00) 
## Version: 
## Last-Updated: Jun 26 2024 (07:24) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(prodlim)
library(survival)
library(randomForestSRC)
library(data.table)
library(rms)
library(riskRegression)
context("riskRegression")
# {{{ survival outcome, pec, Brier Score, external prediction
test_that("survival outcome, Brier Score, external prediction",{
    if (!requireNamespace("pec",quietly=FALSE)){
        message("Package pec not installed. Skip this test.")
        q <- p <- 1
    }else{    
        ## generate simulated data
        set.seed(130971)
        n <- 73
        dat <- sampleData(n,outcome="survival")
        dat[,status:=event]
        setorder(dat,time,-status)
        ## define models
        ## Models <- list("Cox.X1" = coxph(Surv(time,status)~X1,data=dat,y=TRUE),
        Models <- list(
            constant=matrix(rep(0.43,n),ncol=1),
            "runif" = matrix(runif(n),ncol=1),"another.runif" = matrix(runif(n),ncol=1))
        ModelsR <- lapply(Models,function(x)1-x)
        ## training error
        library(prodlim)
        a <- pec::pec(Models,formula = Surv(time,status)~X1+X2,data=dat,times= c(5),exact=FALSE,start=NULL,verbose=TRUE)
        ## compare models
        b <- Score(ModelsR,formula = Surv(time,status)~X1+X2,data=dat,times= c(5),se.fit=FALSE,metrics = "brier")
        q <- b$Brier$score[,Brier]
        p <- as.vector(unlist(a$AppErr))
    }
    expect_equal(ignore_attr=TRUE,q,p)
})

# }}}
# {{{ integrated Brier score

test_that("integrated Brier score",{
    if (!requireNamespace("pec",quietly=TRUE)){
        message("Package pec not installed. Skip this test.")
        p <- q <- 1
    }else{
        set.seed(18)
        trainSurv <- sampleData(100,outcome="survival")
        testSurv <- sampleData(40,outcome="survival")
        library(pec)
        cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
        cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
        ttt = sort(unique(testSurv$time))
        # remove last obs to avoid NA predictions
        ttt = ttt[-length(ttt)]
        suppressWarnings(xs <- Score(list("c1"=cox1,"c2"=cox2),
                                     formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,
                                     se.fit=FALSE,
                                     summary="ibs",
                                     times=ttt))
        library(prodlim)
        xp=pec::pec(list("c1"=cox1,"c2"=cox2),
                    formula=Surv(time,event)~1,data=testSurv,
                    times=ttt)
        a1 <- as.numeric(ibs(xp,times=ttt,models="c1"))
        b1 <- xs$Brier$score[model=="c1"][["IBS"]]
        q <- as.numeric(c(a1,use.names=FALSE))
        p <- c(b1)
    }
    expect_equal(ignore_attr=TRUE,p,q)
})

# }}}
# {{{ survival outcome uncensored

test_that("survival outcome uncensored",{
    if (!requireNamespace("randomForestSRC",quietly=TRUE)){
        message("Package randomForestSRC not installed. Skip this test.")
    }else{
        set.seed(8)
        d <- sampleData(100,outcome="survival")
        d[,event:=rep(1,.N)]
        cx=coxph(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,d,x=TRUE)
        cx1=cph(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,d,x=TRUE,y=TRUE)
        out <- Score(list(Cox=cx,Cph=cx1),data=d,metrics="brier",summary="ibs",contrasts=FALSE,times=1:8,formula=Hist(time,event)~1,se.fit=FALSE)
        expect_equal(ignore_attr=TRUE,100*out$Brier$score[model=="Cox"]$IBS,100*out$Brier$score[model=="Cph"]$IBS,tolerance=0.001)
    }
})

# }}}
# {{{ survival outcome: robustness against order of data set
test_that("survival outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s1 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    a <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="brier")
    A <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="brier",cens.model = "cox")
    setkey(d,X4)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s2 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    setorder(d,time,-event)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s3 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    b <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="brier")
    B <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="brier",cens.model = "cox")
    expect_equal(ignore_attr=TRUE,s1,s2)
    expect_equal(ignore_attr=TRUE,s1,s3)
    expect_equal(ignore_attr=TRUE,a,b)
    expect_equal(ignore_attr=TRUE,A,B)
})
# }}}
# {{{ competing risks outcome: check against pec
test_that("competing risks outcome: check against pec",{
    if (!requireNamespace("pec",quietly=TRUE)){
        message("Package pec not installed. Skip this test.")
    }else{
        set.seed(112)
        d <- sampleData(43,outcome="competing.risks")
        nd <- sampleData(43,outcome="competing.risks")
        f <- FGR(Hist(time,event)~X1+X6,data=d,cause=1)
        a <- pec::pec(list(f),data=nd,times=c(2,5),formula=Hist(time,event)~1,cens.model="marginal",exact=FALSE)
        b <- Score(list(FGR=f),data=nd,formula=Hist(time,event)~1,cens.model="km",se.fit=FALSE,times=c(2,5),metrics="brier")
        expect_equal(ignore_attr=TRUE,a$AppErr$Reference[-1],b$Brier$score[model=="Null model",Brier])
        expect_equal(ignore_attr=TRUE,a$AppErr$FGR[-1],b$Brier$score[model=="FGR",Brier])
    }
})
# }}}
# {{{ competing risks outcome: robustness against order of data set
test_that("competing risks outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="competing.risks")
    f1 <- CSC(Hist(time,event)~X1+X5+X8,data=d)
    f2 <- FGR(Hist(time,event)~X2+X6+X9+X10,data=d,cause=1)
    f3 <- cbind(d$X8,d$X8,d$X8)
    suppressWarnings(s1 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=TRUE,cause=1))
    s1b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=.95,cause=1,metrics="auc")
    setkey(d,X4)
    f3 <- cbind(d$X8,d$X8,d$X8)
    suppressWarnings(s2 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=.95,cause=1))
    s2b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=.95,cause=1,metrics="auc")
    setorder(d,time,-event)
    f3 <- cbind(d$X8,d$X8,d$X8)
    suppressWarnings(s3 <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=.95,cause=1))
    s3b <- Score(list(f1,f2,f3),formula=Hist(time,event)~1,data=d,times=c(3,5,9),conf.int=.95,cause=1,metrics="auc")
    s1$call$conf.int <- .95
    expect_equal(ignore_attr=TRUE,s1,s2)
    expect_equal(ignore_attr=TRUE,s1,s3)
    expect_equal(ignore_attr=TRUE,s1$AUC,s1b$AUC)
    expect_equal(ignore_attr=TRUE,s2$AUC,s2b$AUC)
    expect_equal(ignore_attr=TRUE,s3$AUC,s3b$AUC)
})
# }}}
# {{{ survival outcome: Brier Score pec vs Score
if (requireNamespace("pec",quietly=TRUE)){
    test_that("survival outcome: Brier Score pec vs Score",{
        library(pec)
        set.seed(112)
        d <- sampleData(43,outcome="survival")
        f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
        f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
        p1 <- pec(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,9),exact=FALSE,start=NULL)
        s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,9),conf.int=FALSE,metrics="brier")
        expect_equal(ignore_attr=TRUE,p1$AppErr$coxph,s1$Brier$score[model=="coxph",Brier])
        expect_equal(ignore_attr=TRUE,p1$AppErr$coxph.1,s1$Brier$score[model=="coxph.1",Brier])
        expect_equal(ignore_attr=TRUE,p1$AppErr$Reference,s1$Brier$score[model=="Null model",Brier])
    })
}
# }}}
# {{{ survival outcome: matrix input
test_that("survival outcome: matrix input",{
    set.seed(112)
    dtrain <- sampleData(43,outcome="survival")
    dtest <- sampleData(4,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=dtrain, x = TRUE, y = TRUE)
    f2 <- predictRisk(f1,newdata=dtest,times=c(1:2))
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=dtest,times=c(1:2),conf.int=FALSE,null.model=0L,metrics="brier")
    expect_equal(ignore_attr=TRUE,s1$Brier$score[model=="coxph",Brier],s1$Brier$score[model=="matrix",Brier])
})

# }}}
# {{{ Leave one out bootstrap: Number of models and time points

test_that("Number of models and time points", {
    if (!requireNamespace("pec",quietly=TRUE)){
        message("Package pec not installed. Skip this test.")
    }else{
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
        R1 <- Score(list(a=fit2),data=GBSG2.test,times=1000,B=100,split.method="loob",formula=Surv(time,cens)~1,plots="cali")
        setorder(GBSG2,time,cens)
        ## setorder(GBSG2.test,age)
        GBSG2 <- 7
        r2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(100,500,2000,1000),formula=Surv(time,cens)~1,plots="cali")
        set.seed(11)
        R2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(1000),B=100,split.method="loob",formula=Surv(time,cens)~1,plots="cali")
        ## r1$Calibration$plotframe
        ## r2$Calibration$plotframe[times==1000&model=="a"]
        ## r3 <- pec(list(a=fit2,b=fit1),data=GBSG2.test,exact=FALSE,times=c(1000),formula=Surv(time,cens)~1)
        expect_equal(ignore_attr=TRUE,as.numeric(r1$Brier$score[model=="a"]),
                     as.numeric(r2$Brier$score[model=="a" & times==1000]))
        ## expect_equal(ignore_attr=TRUE,r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    }
})

# }}}
# {{{ Bootstrap cross validation
test_that("Number of models and time points", {
    if (!requireNamespace("pec",quietly=TRUE)){
        message("Package pec not installed. Skip this test.")
    }else{
        library(pec)
        library(survival)
        library(rms)
        library(data.table)
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
        R1 <- Score(list(a=fit2),data=GBSG2.test,seed = 11,times=1000,B=50,split.method="bootcv",formula=Surv(time,cens)~1,plots="cali")
        setorder(GBSG2,time,cens)
        ## setorder(GBSG2.test,age)
        GBSG2 <- 7
        r2 <- Score(list(a=fit2,b=fit1),data=GBSG2.test,times=c(100,500,2000,1000),formula=Surv(time,cens)~1,plots="cali")
        R2 <- Score(list(a=fit2,b=fit1),seed = 11,data=GBSG2.test,times=c(1000),B=50,split.method="bootcv",formula=Surv(time,cens)~1,plots="cali")
        ## r1$Calibration$plotframe
        ## r2$Calibration$plotframe[times==1000&model=="a"]
        ## r3 <- pec(list(a=fit2,b=fit1),data=GBSG2.test,exact=FALSE,times=c(1000),formula=Surv(time,cens)~1)
        expect_equal(ignore_attr=TRUE,as.numeric(r1$Brier$score[model=="a"]),as.numeric(r2$Brier$score[model=="a" & times==1000]))
        expect_equal(ignore_attr=TRUE,as.numeric(R1$Brier$score[model=="a"]),as.numeric(R2$Brier$score[model=="a" & times==1000]))
        ## expect_equal(ignore_attr=TRUE,r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    }
})
# }}}
# {{{ LOOB: Number of models and time points
test_that("LOOB: Number of models and time points", {   
    data(GBSG2,package = "pec")
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
    expect_equal(ignore_attr=TRUE,a$Brier$score[model=="fit1"&times==2000],b$Brier$score[model=="fit1"&times==2000])
    expect_equal(ignore_attr=TRUE,a$Brier$score[model=="fit1"&times==1000],b$Brier$score[model=="fit1"&times==1000])
    expect_equal(ignore_attr=TRUE,a$Brier$score[model=="fit2"&times==2000],b$Brier$score[model=="fit2"&times==2000])
    expect_equal(ignore_attr=TRUE,a$Brier$score[model=="fit2"&times==1000],b$Brier$score[model=="fit2"&times==1000])
    B <- 50
    setorder(GBSG2,time,-cens)
    set.seed(138)
    A <- Score(list(fit1=fit1,fit2=fit2,fit3=fit3),cens.model = "cox",data = GBSG2,times = c(1000,365.25*4),metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
    ## print(A$Brier$score[model=="fit1"&times==365.25*4])
    setorder(GBSG2,time,cens)
    set.seed(138)
    A1 <- Score(list(fit1=fit1),cens.model = "cox",data = GBSG2,times = 365.25*4,metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
    ## print(A1$Brier$score[model=="fit1"&times==365.25*4])
    setorder(GBSG2,age)
    set.seed(138)
    A2 <- Score(list(fit2=fit2,fit1=fit1),cens.model = "cox",data = GBSG2,times = c(365.25*4,300,1000),metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
    ## print(A2$Brier$score[model=="fit1"&times==365.25*4])
    expect_equal(ignore_attr=TRUE,A$Brier$score[model=="fit1"&times==365.25*4],
                 A2$Brier$score[model=="fit1"&times==365.25*4])
    expect_equal(ignore_attr=TRUE,A$Brier$score[model=="fit1"&times==1000],
                 A2$Brier$score[model=="fit1"&times==1000])
    expect_equal(ignore_attr=TRUE,A$Brier$score[model=="fit2"&times==1000],
                 A2$Brier$score[model=="fit2"&times==1000])
    expect_equal(ignore_attr=TRUE,A1$Brier$score[model=="fit1"&times==365.25*4],
                 A2$Brier$score[model=="fit1"&times==365.25*4])
    expect_equal(ignore_attr=TRUE,A1$Brier$score[model=="fit1"&times==365.25*4],
                 A$Brier$score[model=="fit1"&times==365.25*4])
})
# }}}
# {{{ vcov AUC
test_that("vcov AUC",{
    set.seed(112)
    d <- sampleData(43,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    test <- Score(list(f1,f2),keep="vcov",formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
    expect_equal(ignore_attr=TRUE,dim(test$AUC$vcov),c(2,2))
    ## survival
    set.seed(112)
    d <- sampleData(112,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x=TRUE,y=TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x=TRUE,y=TRUE)
    test <- Score(list(a=f1,f2),times=c(5,7),keep="vcov",formula=Surv(time,event)~1,data=d,conf.int=TRUE,metrics=c("brier","auc"))
    expect_equal(ignore_attr=TRUE,dim(test$AUC$vcov),c(4,4))
})
# }}}
# {{{ LOOB binary
## GIVES WARNING 
## already exporting variable(s): data, split.method, Weights, N, trainseeds
test_that("loob binary",{
    learndat=sampleData(200,outcome="binary")
    lr1a = glm(Y~X6,data=learndat,family=binomial)
    lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
    ## leave-one-out bootstrap
    ## in series
    set.seed(5)
    system.time(loob.se0 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=TRUE))
    ## in multicore 3 cpu's 
    set.seed(5)
    system.time(loob.se0a <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,parallel="multicore",ncpus=3,se.fit=TRUE))
    ## snow 3 cpu's 
    library(parallel)
    mycl <- parallel::makePSOCKcluster(rep("localhost", 3))
    set.seed(5)
    system.time(loob.se0b <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,parallel="snow",cl=mycl,se.fit=TRUE))
    try(stopCluster(mycl),silent=TRUE)
    set.seed(5)
    loob.se1 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=FALSE)
    expect_equal(ignore_attr=TRUE,loob.se0$AUC$contrasts$delta,loob.se0b$AUC$contrasts$delta)
    expect_equal(ignore_attr=TRUE,loob.se0$AUC$contrasts$delta,loob.se0a$AUC$contrasts$delta)
    expect_equal(ignore_attr=TRUE,loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
    expect_equal(ignore_attr=TRUE,loob.se0$Brier$contrasts$delta,loob.se1$Brier$contrasts$delta)
})
# }}}
# {{{ bootcv binary
test_that("bootcv binary (multi.state.test)",{
    learndat=sampleData(200,outcome="binary")
    lr1a = glm(Y~X6,data=learndat,family=binomial)
    lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
    ## leave-one-out bootstrap
    set.seed(5)
    bootcv.se0 <- Score(list("LR1"=lr1a,"LR2"=lr2a),conservative=TRUE,formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=FALSE,metric="brier")
    set.seed(5)
    bootcv.se1 <- Score(list("LR1"=lr1a,"LR2"=lr2a),conservative=TRUE,formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=FALSE,metric="brier")
    set.seed(5)
    bootcv.se2 <- Score(list("LR1"=lr1a,"LR2"=lr2a),conservative=TRUE,formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=FALSE,multi.split.test=TRUE,metric="brier")
    set.seed(5)
    bootcv.se3 <- Score(list("LR1"=lr1a,"LR2"=lr2a),conservative=TRUE,formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=TRUE,multi.split.test=TRUE,metric="brier")
    bootcv <- list(bootcv.se0,bootcv.se1,bootcv.se2,bootcv.se3)
    ## delta
    for (i in 1:4)
        for (j in 2:4) 
            for (m in c("Brier"))
                expect_equal(ignore_attr=TRUE,bootcv[[i]][[m]]$contrasts$delta,bootcv[[j]][[m]]$contrasts$delta)
    ## lower, upper
    expect_equal(ignore_attr=TRUE,bootcv[[2]][["Brier"]]$contrasts[,.(lower,upper)],bootcv[[4]][["Brier"]]$contrasts[,.(lower,upper)])
    ## p-value (does not work yet)
    ## for (m in c("AUC","Brier")){
    ## expect_equal(ignore_attr=TRUE,bootcv[[3]][[m]]$contrasts[,.(p)],bootcv[[4]][[m]]$contrasts[,.(p)])
    ## }
})
# }}}
# {{{ "Brier score"
if (requireNamespace("pec",quietly=TRUE)){
    test_that("Brier score",{
        library(pec)
        library(riskRegression)
        library(survival)
        data(Melanoma)
        fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
        ## predictRisk(fit.lrr,times=c(1,10,100,1000),newdata=Melanoma)
        fit.arr2 <- ARR(Hist(time,status)~thick+age,data=Melanoma,cause=1)
        fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
        ## fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
        system.time(old <- pec(list(
                        ARR=fit.arr2a,
                        ## ARR.power=fit.arr2b,
                        LRR=fit.lrr),
                        data=Melanoma,
                        formula=Hist(time,status)~1,
                        times = c(500,1000,2000,3000),
                        exact = FALSE,
                        start = NULL,
                        cause=1, B=10,split.method="none"))
        ## predictRisk(fit.arr2a,newdata=Melanoma[1:10,],times=0)
        system.time(new <- Score(list(
                        ARR=fit.arr2a,
                        ## ARR.power=fit.arr2b,
                        LRR=fit.lrr),
                        data=Melanoma,se.fit = FALSE,
                        ## times=c(0,sort(unique(Melanoma$time))),
                        metrics="brier",plots=NULL,summary=NULL,
                        formula=Hist(time,status)~1,
                        times = c(500,1000,2000,3000),
                        cause=1, B=10,split.method="none"))
        nix <- lapply(1:3,function(m){
            expect_equal(ignore_attr=TRUE,
                         new$Brier$score[model==names(new$models)[m]][["Brier"]],
                         old$AppErr[[names(old$AppErr)[[m]]]])})
    })
}
# }}}
# {{{ cv-k
test_that("loob survival",{
    set.seed(8)
    learndat=sampleData(188,outcome="survival")
    cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
    cox2a = coxph(Surv(time,event)~X1+X9,data=learndat,x=TRUE,y=TRUE)
    ## leave-one-out bootstrap
    x <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="loob",B=100,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    y <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="loob",B=10,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    z <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="loob",M = .632*nrow(learndat),B=100,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    ## 10-fold, 7-fold, 5-fold, 2-fold
    a <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="cv10",B=1,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    b <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="cv7",B=2,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    c <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="cv2",B=20,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    d <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="cv5",B=2,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    ## bootcv
    e <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="bootcv",B=100,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
})
# }}}
# {{{ sens, spec, ppv, npv against timeROC
testthat("sens, spec, ppv, npv against timeROC",{
    library(timeROC)
    data(pbc)
    pbc<-pbc[!is.na(pbc$trt),] # select only randomised subjects
    pbc$status<-as.numeric(pbc$status==2) # create event indicator: 1 for death, 0 for censored     
    # Se, Sp, PPV and NPV computation for serum bilirunbin at threshold c=0.9(mg/dl) 
    res.SeSpPPVNPV.bili <- SeSpPPVNPV(cutpoint=0.9,
                                      T=pbc$time,
                                      delta=pbc$status,marker=pbc$bili,
                                      cause=1,weighting="marginal",
                                      times=1200,
                                      iid=TRUE)
    x <- Score(list(pbc$bili),formula=Surv(time,status)~1,data=pbc,times=1200,metrics="auc", cutpoints=0.9)
    ## Check that they give the same results
    a = as.numeric(x$AUC$cutpoints[,c(4,6,8,10)])
    b = as.numeric(sapply(res.SeSpPPVNPV.bili[c(1,2,3,4)],"[",2))
    expect_equal(a,b)
    se.a = as.numeric(x$AUC$cutpoints[,c(5,7,9,11)])
    se.b = as.numeric(sapply(res.SeSpPPVNPV.bili$inference[c(7,9,10,12)],"[",2))
    expect_equal(se.a,se.b)    

    ##-------------With competing risks-------------------
    ## Also test computation time
    data(Paquid)
    system.time(res.SeSpPPVNPV.DSST <- SeSpPPVNPV(cutpoint=22,
                                                  T=Paquid$time,
                                                  delta=Paquid$status,marker=Paquid$DSST,
                                                  cause=1,weighting="marginal",
                                                  times=5,iid=TRUE))
    library(riskRegression)
    system.time(x<-Score(list(Paquid$DSST),formula=Hist(time,status)~1,data=Paquid,times=c(5),metrics="auc", cutpoints=22))
    ## Check that they give the same results
    a = as.numeric(x$AUC$cutpoints[,c(4,6,8,10)])
    b = as.numeric(sapply(res.SeSpPPVNPV.DSST[c(1,3,4,6)],"[",2))
    expect_equal(a,b)
    a.se = as.numeric(x$AUC$cutpoints[,c(5,7,9,11)])
    b.se = as.numeric(sapply(res.SeSpPPVNPV.DSST$inference[c(7,9,10,12)],"[",2))
    expect_equal(a.se,b.se,tolerance = 0.0001)
    ## Test Cox censoring
    ## x<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
    ## y<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22,censoring.save.memory = TRUE)
    ## rbind(x$AUC$cutpoints,
    ## y$AUC$cutpoints)
})
# }}}
# {{{ cutpoints with competing risks
testthat("cutpoints with competing risks",{
    ## Also test computation time
    data(Paquid,package = "timeROC")
    system.time(res.SeSpPPVNPV.DSST <- SeSpPPVNPV(cutpoint=22,
                                                  T=Paquid$time,
                                                  delta=Paquid$status,marker=Paquid$DSST,
                                                  cause=1,weighting="marginal",
                                                  times=5,iid=TRUE))

library(riskRegression)
system.time(x<-Score(list(Paquid$DSST),formula=Hist(time,status)~1,data=Paquid,times=c(5),metrics="auc", cutpoints=22))
## Check that they give the same results
x$AUC$res.cut[,c(4,6,8,10)]
res.SeSpPPVNPV.DSST[c(1,3,4,6)]
x$AUC$res.cut[,c(5,7,9,11)]
res.SeSpPPVNPV.DSST$inference[c(7,9,10,12)]

## Test Cox censoring
x<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
x$AUC$res.cut
y<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22,censoring.save.memory = TRUE)
y$AUC$res.cut
## They are a bit different.
}
# }}}
######################################################################
### test-Score_survival_outcome.R ends here
