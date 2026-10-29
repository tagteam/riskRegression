library(testthat)
library(data.table)
context("Prediction error")
# {{{ "AUC censored data: order"
test_that("AUC censored data: order",{
    library(riskRegression)
    library(survival)
    library(prodlim)
    data(Melanoma)
    setDT(Melanoma)
    setkey(Melanoma,age)
    fit <- coxph(Surv(time,status!=0)~invasion+epicel+logthick,data=Melanoma,x=TRUE)
    suppressWarnings(a <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="marginal",metric="Auc"))
    suppressWarnings(A <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="cox",metric="Auc"))
    setkey(Melanoma,logthick)
    suppressWarnings(b <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="marginal",metric="Auc"))
    suppressWarnings(B <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="cox",metric="Auc"))
    a$call <- b$call <- A$call <- B$call <- NULL
    expect_equal(ignore_attr=TRUE,A,B)
    expect_equal(ignore_attr=TRUE,a,b)
})
# }}}
# {{{ "AUC: standard error for n > 46340 (integer overflow of n*n)"
## 4 copies of 12500 subjects give n=50000 > sqrt(2^31-1).
## The AUC is unchanged by replication and its standard error should halve.
test_that("AUC: standard error for n > 46340, cens.model km",{
    library(riskRegression)
    library(survival)
    library(prodlim)
    set.seed(8)
    n <- 12500
    X <- rnorm(n)
    T <- rexp(n,0.1*exp(0.7*X))
    C <- rexp(n,0.03)
    d <- data.table(time=pmin(T,C),status=as.integer(T<=C),X=X)
    d4 <- d[rep(1:n,4)]
    for (conservative in c(TRUE,FALSE)){
        a <- Score(list(X=d$X),data=d,Hist(time,status)~X,times=5,metrics="auc",null.model=FALSE,cens.model="km",conservative=conservative)$AUC$score
        a4 <- Score(list(X=d4$X),data=d4,Hist(time,status)~X,times=5,metrics="auc",null.model=FALSE,cens.model="km",conservative=conservative)$AUC$score
        expect_equal(a4$AUC,a$AUC)
        expect_equal(2*a4$se,a$se,tolerance=1e-3)
    }
})
test_that("AUC: standard error for n > 46340, cens.model cox",{
    skip_on_cran()
    library(riskRegression)
    library(survival)
    library(prodlim)
    set.seed(8)
    n <- 12500
    X <- rnorm(n)
    T <- rexp(n,0.1*exp(0.7*X))
    C <- rexp(n,0.03)
    d <- data.table(time=pmin(T,C),status=as.integer(T<=C),X=X)
    d4 <- d[rep(1:n,4)]
    a <- Score(list(X=d$X),data=d,Hist(time,status)~X,times=5,metrics="auc",null.model=FALSE,cens.model="cox",conservative=FALSE,censoring.save.memory=TRUE)$AUC$score
    a4 <- Score(list(X=d4$X),data=d4,Hist(time,status)~X,times=5,metrics="auc",null.model=FALSE,cens.model="cox",conservative=FALSE,censoring.save.memory=TRUE)$AUC$score
    expect_false(is.na(a4$se))
    ## not exact: ties from replication change the Cox model for the censoring distribution
    expect_equal(2*a4$se,a$se,tolerance=0.05)
})
# }}}

