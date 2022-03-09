### slowtest-score.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Dec  6 2020 (09:25) 
## Version: 
## Last-Updated: Mar  9 2022 (08:27) 
##           By: Thomas Alexander Gerds
##     Update #: 7
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
library(riskRegression)
library(data.table)

# {{{ "survival outcome: robustness against order of data set"
test_that("survival outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
    f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s1 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=TRUE,metrics="auc")
    setkey(d,X4)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s2 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    setorder(d,time,-event)
    f3 <- cbind(d$X8,d$X8,d$X8)
    s3 <- Score(list(f1,f2,f3),formula=Surv(time,event)~1,data=d,times=c(3,5,10),conf.int=.95,metrics="auc")
    s1$call$conf.int <- .95
    expect_equal(s1,s2)
    expect_equal(s1,s3)
})
# }}}

# {{{ "competing risks outcome: check against pec"
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
        expect_equal(a$AppErr$Reference[-1],b$Brier$score[model=="Null model",Brier])
        expect_equal(a$AppErr$FGR[-1],b$Brier$score[model=="FGR",Brier])
    }
})

# }}}
# {{{ "competing risks outcome: robustness against order of data set"
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
    expect_equal(s1,s2)
    expect_equal(s1,s3)
    expect_equal(s1$AUC,s1b$AUC)
    expect_equal(s2$AUC,s2b$AUC)
    expect_equal(s3$AUC,s3b$AUC)
})
# }}}

# {{{ "survival outcome: Brier Score pec vs Score"
if (requireNamespace("pec",quietly=TRUE)){
    test_that("survival outcome: Brier Score pec vs Score",{
        library(pec)
        set.seed(112)
        d <- sampleData(43,outcome="survival")
        f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d, x = TRUE, y = TRUE)
        f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d, x = TRUE, y = TRUE)
        p1 <- pec(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,9),exact=FALSE,start=NULL)
        s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,9),conf.int=FALSE,metrics="brier")
        expect_equal(p1$AppErr$coxph,s1$Brier$score[model=="coxph",Brier])
        expect_equal(p1$AppErr$coxph.1,s1$Brier$score[model=="coxph.1",Brier])
        expect_equal(p1$AppErr$Reference,s1$Brier$score[model=="Null model",Brier])
    })
}
# }}}

# {{{ "survival outcome: matrix input"
test_that("survival outcome: matrix input",{
    set.seed(112)
    dtrain <- sampleData(43,outcome="survival")
    dtest <- sampleData(4,outcome="survival")
    f1 <- coxph(Surv(time,event)~X1+X5+X8,data=dtrain, x = TRUE, y = TRUE)
    f2 <- predictRisk(f1,newdata=dtest,times=c(1:2))
    s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=dtest,times=c(1:2),conf.int=FALSE,null.model=0L,metrics="brier")
    expect_equal(s1$Brier$score[model=="coxph",Brier],s1$Brier$score[model=="matrix",Brier])
})

# }}}
# {{{ "Leave one out bootstrap: Number of models and time points"
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
        expect_equal(as.numeric(r1$Brier$score[model=="a"]),
                     as.numeric(r2$Brier$score[model=="a" & times==1000]))
        ## expect_equal(r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    }
})

# {{{ "Bootstrap cross validation
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
        expect_equal(as.numeric(r1$Brier$score[model=="a"]),as.numeric(r2$Brier$score[model=="a" & times==1000]))
        ## expect_equal(r1$AUC$score[model=="a"],r2$AUC$score[model=="a" & times==1000])
    }
})
# }}}

# {{{ "LOOB: Number of models and time points"
test_that("LOOB: Number of models and time points", {   
    library(testthat)
    library(survival)
    library(rms)
    library(riskRegression)
    library(data.table)
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
    ## print(A$Brier$score[model=="fit1"&times==365.25*4])
    setorder(GBSG2,time,cens)
    set.seed(138)
    A1 <- Score(list(fit1=fit1),cens.model = "cox",data = GBSG2,times = 365.25*4,metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
    ## print(A1$Brier$score[model=="fit1"&times==365.25*4])
    setorder(GBSG2,age)
    set.seed(138)
    A2 <- Score(list(fit2=fit2,fit1=fit1),cens.model = "cox",data = GBSG2,times = c(365.25*4,300,1000),metric = "brier",null.model = FALSE,contrast = FALSE,B=B,split.method="looboot",formula = Surv(time, cens)~horTh)
    ## print(A2$Brier$score[model=="fit1"&times==365.25*4])
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



# {{{ "vcov AUC"
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
# }}}


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
    expect_equal(loob.se0$AUC$contrasts$delta,loob.se0b$AUC$contrasts$delta)
    expect_equal(loob.se0$AUC$contrasts$delta,loob.se0a$AUC$contrasts$delta)
    expect_equal(loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
    expect_equal(loob.se0$Brier$contrasts$delta,loob.se1$Brier$contrasts$delta)
})

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
                expect_equal(bootcv[[i]][[m]]$contrasts$delta,bootcv[[j]][[m]]$contrasts$delta)
    ## lower, upper
    expect_equal(bootcv[[2]][["Brier"]]$contrasts[,.(lower,upper)],bootcv[[4]][["Brier"]]$contrasts[,.(lower,upper)])
    ## p-value (does not work yet)
    ## for (m in c("AUC","Brier")){
    ## expect_equal(bootcv[[3]][[m]]$contrasts[,.(p)],bootcv[[4]][[m]]$contrasts[,.(p)])
    ## }
})

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
        fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
        system.time(old <- pec(list(ARR=fit.arr2a,ARR.power=fit.arr2b,LRR=fit.lrr),
                               data=Melanoma,
                               formula=Hist(time,status)~1,
                               cause=1, B=10,split.method="none"))
        ## predictRisk(fit.arr2a,newdata=Melanoma[1:10,],times=0)
        system.time(new <- Score(list(ARR=fit.arr2a,ARR.power=fit.arr2b,LRR=fit.lrr),
                                 data=Melanoma,conf.int=0,
                                 times=c(0,sort(unique(Melanoma$time))),
                                 metrics="brier",plots=NULL,summary=NULL,
                                 formula=Hist(time,status)~1,
                                 cause=1, B=10,split.method="none"))
        nix <- lapply(1:4,function(m){
            expect_equal(new$Brier$score[model==names(new$models)[m]][["Brier"]],
                         old$AppErr[[names(old$AppErr)[[m]]]])})
    })
}
# }}}


######################################################################
### slowtest-score.R ends here
