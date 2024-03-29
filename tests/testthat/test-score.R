### test-Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (14:30) 
## Version: 
## last-updated: Sep 17 2022 (07:00) 
##           By: Thomas Alexander Gerds
##     Update #: 172
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
context("riskRegression")

# {{{ Missing values
test_that("Missing values in data", {
    d <- data.frame(time=c(1,2,3,NA,5,6),event=c(1,NA,1,0,NA,0),X=c(1,3,1,NA,9,-8))
    expect_error(Score(list(d$X),data=d,times=3,formula=Hist(time,event)~1,metrics="auc"))
})
# }}}
# {{{ "R squared/IPA"
test_that("R squared/IPA", { 
    set.seed(112)
    d <- sampleData(43,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    r1 <- rsquared(f1,newdata=d)
    r2 <- IPA(f2,newdata=d)
    full <- Score(list(f1=f1,f2=f2),formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
    expect_equal(ignore_attr=TRUE,r1$IPA.drop[1],full$Brier$score[model=="f1",IPA])
    expect_equal(ignore_attr=TRUE,r2$IPA[2],full$Brier$score[model=="f2",IPA])
})
# }}}

# {{{ "binary outcome: robustness against order of data set"
test_that("binary outcome: robustness against order of data set",{
    set.seed(112)
    d <- sampleData(43,outcome="binary")
    f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
    f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
    f3 <- d$X8
    s1 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=TRUE,metrics="auc")
    s1b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    setkey(d,X4)
    f3 <- d$X8
    s2 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    s2b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    setorder(d,Y)
    f3 <- d$X8
    s3 <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    s3b <- Score(list(f1,f2,f3),formula=Y~1,data=d,conf.int=.95,metrics="auc")
    ## lapply(names(s1),function(n){print(n);expect_equal(ignore_attr=TRUE,s1[[n]],s3[[n]])})
    s1$call$conf.int <- .95
    expect_equal(ignore_attr=TRUE,s1,s2)
    expect_equal(ignore_attr=TRUE,s1,s3)
    expect_equal(ignore_attr=TRUE,s1$AUC,s1b$AUC)
    expect_equal(ignore_attr=TRUE,s2$AUC,s2b$AUC)
    expect_equal(ignore_attr=TRUE,s3$AUC,s3b$AUC)
})
# }}}
# {{{ "survival outcome, Brier Score, external prediction"
test_that("survival outcome,Brier Score, external prediction",{
    if (!requireNamespace("pec",quietly=FALSE)){
        message("Package pec not installed. Skip this test.")
        q <- p <- 1
    }else{    
        ## generate simulated data
        set.seed(130971)
        n <- 100
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
        b <- Score(ModelsR,formula = Surv(time,status)~X1+X2,data=dat,times= c(5),se.fit=FALSE)
        q <- b$Brier$score[,Brier]
        p <- as.vector(unlist(a$AppErr))
    }
    expect_equal(ignore_attr=TRUE,q,p)
})

# }}}
# {{{integrated Brier score
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
        ## cbind(a1,b1)
        q <- as.numeric(c(a1,use.names=FALSE))
        p <- c(b1)
    }
    expect_equal(ignore_attr=TRUE,p,q)
})
# }}}

# {{{ "survival outcome uncensored"
test_that("survival outcome uncensored",{
    if (!requireNamespace("randomForestSRC",quietly=TRUE)){
        message("Package randomForestSRC not installed. Skip this test.")
    }else{
        library(survival)
        library(data.table)
        library(rms)
        library(riskRegression)
        library(prodlim)
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
    expect_equal(ignore_attr=TRUE,s1,S1)
    expect_equal(ignore_attr=TRUE,s2,S2)
    expect_equal(ignore_attr=TRUE,s3,S3)
    expect_equal(ignore_attr=TRUE,s1$Brier,s2$Brier)
    expect_equal(ignore_attr=TRUE,s1$Brier$score$Brier,s3$Brier$score$Brier)
    expect_equal(ignore_attr=TRUE,s2$Brier$score$Brier,s3$Brier$score$Brier)
    expect_equal(ignore_attr=TRUE,S1$Brier,S2$Brier)
    expect_equal(ignore_attr=TRUE,S1$Brier$score$Brier,S3$Brier$score$Brier)
    expect_equal(ignore_attr=TRUE,S2$Brier$score$Brier,S3$Brier$score$Brier)
})
# }}}
# {{{ "binary outcome: AUC"
test_that("binary outcome: AUC", {
    if (!requireNamespace("pROC",quietly=TRUE)){
        message("Package pROC not installed. Skip this test. predictCSC.")
    }else{
        set.seed(17)
        y <- rbinom(100, 1, .5)
        x1 <- rnorm(100) + 1.5 * y
        x2 <- rnorm(100) + .5 * y
        x3 <- rnorm(100) + 2.5 * y
        x <- data.frame(x1,x2,x3)
        y <- as.factor(y)
        r1 <- pROC::roc(y~x1)
        r2 <- pROC::roc(y~x2)
        r3 <- pROC::roc(y~x3)
        procres <- pROC::roc.test(r1,r2)
        d <- data.frame(x1,x2,x3,y)
        ## Source(riskRegression)
        scoreres <- Score(list(X1=~x1,X2=~x2,X3=~x3),formula=y~1,data=d,null.model=FALSE,cause="1",metrics="auc")
        ## Roc(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d)
        scoreres <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,cause="1")
        ## to avoid side effects of data.table features we check the following 
        scoreres1 <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,metrics="auc",cause="1")
        scoreres1a <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,null.model=FALSE,metrics="auc",se.fit=0L,cause="1")
        expect_equal(ignore_attr=TRUE,scoreres$AUC,scoreres1$AUC)
        ## daim.auc <- daimres$AUC[,c("AUC","SD(DeLong)")]
        score.auc <- as.data.frame(scoreres$AUC$score[,c("AUC","se"),with=FALSE])
        ## rownames(score.auc) <- rownames(daim.auc)
        ## colnames(score.auc) <- colnames(daim.auc)
        ## expect_equal(ignore_attr=TRUE,daim.auc,score.auc)
        expect_equal(ignore_attr=TRUE,scoreres$AUC$score[["AUC"]],c(r1$auc,r2$auc,r3$auc))
        score.diff <- scoreres$AUC$contrasts[,c("delta.AUC","se","lower","upper","p"),with=FALSE]
        ## daim.diff <- daimres$difference
        ## expect_equal(ignore_attr=TRUE,daim.diff$"AUC Difference",-score.diff$delta.AUC)
        ## expect_equal(ignore_attr=TRUE,daim.diff$"CI(lower)",-score.diff$upper)
        ## expect_equal(ignore_attr=TRUE,daim.diff$"CI(upper)",-score.diff$lower)
        ## expect_equal(ignore_attr=TRUE,daim.diff$"P.Value",score.diff$p)
    }
})
# }}}

## library(survival)
## library(riskRegression)
## library(rms)
## data(pbc)
## pbc <- na.omit(pbc)
## pbc$time=pbc$time+rnorm(nrow(pbc),sd=.1)
## a <- cph(Surv(time,status!=0)~age+edema+sex+log(bili),data=pbc,surv=TRUE,y=1,x=1)
## b <- cph(Surv(time,status!=0)~age+edema+sex+log(bili)+log(protime)+log(albumin),data=pbc,surv=TRUE,y=1,x=1)
## set.seed(17)
## x <- Score(list(a,b),data=pbc,formula=Surv(time,status!=1)~1,cause=1,times=c(1000),metrics=c("auc"))
## r <- pec(list(a,b),data=pbc,start=NULL,Surv(time,status!=1)~1,times=c(100,500,1000),exact=FALSE)
## u <- with(pbc,timeROC(T=time,delta=status!=0,marker=1-predictSurvProb(a,times=1500,newdata=pbc),cause=1,times=1500,iid=TRUE))
## u2 <- with(pbc,timeROC(T=time,delta=status!=0,marker=1-predictSurvProb(b,times=1500,newdata=pbc),cause=1,times=c(1500)))
## v <- Score(list(a,b),data=pbc,formula=Surv(time,status!=0)~1,times=c(500,1500),metrics=c("AUC"))

#----------------------------------------------------------------------
### test-Score.R ends here
