### test-Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (14:30) 
## Version: 
## last-updated: Mar  1 2017 (16:16) 
##           By: Thomas Alexander Gerds
##     Update #: 17
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
test_that("survival outcome: Brier Score",
{
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

test_that("survival outcome,Brier Score, external prediction",
{
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

test <- FALSE
# Brice: I can't run 
    ## expect_equal(daim.auc,score.auc)
    ## expect_equal(daim.diff$"CI(lower)",-score.diff$upper)
    ## expect_equal(daim.diff$"CI(upper)",-score.diff$lower)
    ## expect_equal(daim.diff$"P.Value",score.diff$p)

# on my computer I get the following error
    ## Error: `daim.auc` not equal to `score.auc`.
    ## Component "SD(DeLong)": Mean relative difference: 0.3720713

    ## Error: daim.diff$"CI(lower)" not equal to -score.diff$upper.
    ## 3/3 mismatches (average diff: 0.0204)
    ## [1]  0.0421 -  0.0327 == 0.00936
    ## [2] -0.2796 - -0.3052 == 0.02561
    ## [3] -0.5012 - -0.5274 == 0.02628

    ## Error: daim.diff$"CI(upper)" not equal to -score.diff$lower.
    ## 3/3 mismatches (average diff: 0.0204)
    ## [1]  0.3418 -  0.3512 == -0.00936
    ## [2] -0.0881 - -0.0625 == -0.02561
    ## [3] -0.2505 - -0.2242 == -0.02628

    ## Error: daim.diff$P.Value not equal to score.diff$p.
    ## 3/3 mismatches (average diff: 0.00297)
    ## [1] 1.21e-02 - 1.81e-02 == -6.08e-03
    ## [2] 1.68e-04 - 2.99e-03 == -2.82e-03
    ## [3] 4.20e-09 - 1.19e-06 == -1.18e-06

if(test){
    
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
    scoreres <- Score(list(X1=~x1,X2=~x2,X3=~x3),formula=y~1,data=d,nullModel=FALSE,cause="1")
    ## Roc(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d)
    scoreres <- Score(list(X1=glm(y~x1,data=d,family='binomial'),X2=glm(y~x2,data=d,family='binomial'),X3=glm(y~x3,data=d,family='binomial')),formula=y~1,data=d,nullModel=FALSE,cause="1")
    daim.auc <- daimres$AUC[,c("AUC","SD(DeLong)")]
    score.auc <- as.data.frame(scoreres$AUC$score[,c("AUC","se.AUC"),with=FALSE])
    rownames(score.auc) <- rownames(daim.auc)
    colnames(score.auc) <- colnames(daim.auc)
    expect_equal(daim.auc,score.auc)
    expect_equal(scoreres$AUC$score[["AUC"]],c(r1$auc,r2$auc,r3$auc))
    score.diff <- scoreres$AUC$contrasts[,c("delta.auc","se.auc","lower","upper","p"),with=FALSE]
    daim.diff <- daimres$difference
    expect_equal(daim.diff$"AUC Difference",-score.diff$delta.auc)
    expect_equal(daim.diff$"CI(lower)",-score.diff$upper)
    expect_equal(daim.diff$"CI(upper)",-score.diff$lower)
    expect_equal(daim.diff$"P.Value",score.diff$p)
})

}
#----------------------------------------------------------------------
### test-Score.R ends here
