library(riskRegression)
library(survival)
library(testthat)
testthat("loob binary",{
    learndat=sampleData(200,outcome="binary")
    lr1a = glm(Y~X6,data=learndat,family=binomial)
    lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
    ## leave-one-out bootstrap
    set.seed(5)
    loob.se0 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=FALSE)
    set.seed(5)
    loob.se1 <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=TRUE)
    expect_equal(loob.se0$AUC$contrasts$delta,loob.se1$AUC$contrasts$delta)
})
## bootstrap cross-validation
x1a=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100,se.fit=FALSE)
x1

x1a=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=10,se.fit=TRUE)
x1a=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="brier",B=10,se.fit=TRUE)
x1b=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="auc",B=10,se.fit=TRUE)
x1c=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="auc",B=10,se.fit=TRUE,multi.split.test=TRUE)

Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="brier",B=10,se.fit=FALSE,multi.split.test=TRUE)
Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="brier",B=10,se.fit=TRUE,multi.split.test=TRUE)
x1d <- Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",metric="brier",B=10,se.fit=TRUE,multi.split.test=FALSE)

x1a=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100,se.fit=TRUE)
x1=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=TRUE)
