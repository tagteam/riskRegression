### test-Score_binary_outcome.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 17 2024 (11:59) 
## Version: 
## Last-Updated: Jul  2 2024 (11:51) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(riskRegression)
library(data.table)
context("binary outcome")
set.seed(112)
d <- sampleData(43,outcome="binary")
f1 <- glm(Y~X1+X5+X8,data=d, family="binomial")
f2 <- glm(Y~X2+X6+X9+X10,data=d, family="binomial")
f3 <- d$X8
# {{{ Missing values in data
test_that("Missing values in data", {
    d <- data.frame(time=c(1,2,3,NA,5,6),event=c(1,NA,1,0,NA,0),X=c(1,3,1,NA,9,-8))
    expect_error(Score(list(d$X),data=d,times=3,formula=Hist(time,event)~1,metrics="auc"))
})
# }}}
# {{{ IPA
test_that("R squared/IPA", { 
    r1 <- rsquared(f1,newdata=d)
    r2 <- IPA(f2,newdata=d)
    full <- Score(list(f1=f1,f2=f2),formula=Y~1,data=d,conf.int=TRUE,summary=c("RR"),plots="ROC")
    expect_equal(ignore_attr=TRUE,r1$IPA.drop[1],full$Brier$score[model=="f1",IPA])
    expect_equal(ignore_attr=TRUE,r2$IPA[2],full$Brier$score[model=="f2",IPA])
})
# }}}
# {{{ robustness against order of data
test_that("robustness against order of data set. metrics = auc",{
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
test_that("binary outcome: robustness against order of data. metrics: Brier and auc + Brier",{
    s1 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics="brier",cause="1")
    s2 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics=c("auc","brier"),cause="1")
    s3 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics=c("auc","brier"),se.fit=FALSE,cause="1")
    setkey(d,Y)
    S1 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics="brier",cause="1")
    S2 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics=c("auc","brier"),cause="1")
    S3 <- Score(list(X6=glm(Y~X6,data=d,family='binomial'),X9=glm(Y~X9,data=d,family='binomial'),X10=glm(Y~X10,data=d,family='binomial')),formula=Y~1,data=d,null.model=FALSE,metrics=c("auc","brier"),se.fit=FALSE,cause="1")
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
# {{{ binary outcome: AUC comparison with pROC
test_that("binary outcome: AUC comparison with pROC", {
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
        score.auc <- as.data.frame(scoreres$AUC$score[,c("AUC","se"),with=FALSE])
        expect_equal(ignore_attr=TRUE,scoreres$AUC$score[["AUC"]],c(r1$auc,r2$auc,r3$auc))
    }
})
# }}}
# {{{ print and summary functionality
test_that("print and summary functionality without null.model", {
    x <- Score(list(d$X8),formula=Y~1,data=d,conf.int=TRUE,metrics="AUC",null.model = FALSE)
    expect_output(print(x))
    expect_output(summary(x))
})
# }}}

# {{{ cutpoints, sens, spec, PPV, NPV
## testthat("cutpoints, sens, spec, PPV, NPV",{
    ## data(pbc,package = "survival")
    ## setDT(pbc)
    ## dd = pbc[,.(Y = 1*(time <1200),bili = bili,X = 1*(bili>0.91))]
    ## Y = dd$Y
    ## X = dd$X
    ## x <- Score(list(dd$bili,dd$X),formula=Y~1,data=dd,metrics="auc",cutpoints = 0.91,breaks = rev(sort(unique(dd$bili))))
    ## u = table(highbili = X,event = Y)
    ## U = as.numeric(u)
    ## names(U) = c("lowbili/noevent","highbili/noevent","lowbili/event","highbili/event")
    ## a = data.table(Sens = U["highbili/event"]/(U["lowbili/event"]+U["highbili/event"]),
          ## Spec = U["lowbili/noevent"]/(U["lowbili/noevent"]+U["highbili/noevent"]),
          ## PPV = U["highbili/event"]/(U["highbili/event"]+U["highbili/noevent"]),
          ## NPV = U["lowbili/noevent"]/(U["lowbili/event"]+U["lowbili/noevent"]))
    ## b = x$AUC$cutpoints[,.(TPR,FPR,PPV,NPV)]
    ## b
    ## a
## })
# }}}

######################################################################
### test-Score_binary_outcome.R ends here
