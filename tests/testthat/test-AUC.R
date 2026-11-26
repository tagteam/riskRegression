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

# {{{ "AUC binary: cutpoints with ties"
## Predicted positive is risk >= cutpoint. The standard errors should use the
## same rule as the estimates, also when risks are tied at the cutpoint.
test_that("AUC binary: standard errors at cutpoints with ties",{
    library(riskRegression)
    set.seed(8)
    d <- data.frame(score = rep(0:4, each = 40))
    d$Y <- rbinom(200, 1, plogis(-2 + 0.8 * d$score))
    x <- Score(list(score = d$score),formula = Y~1,data = d,metrics = "auc",
               cutpoints = 2,null.model = FALSE,verbose = -1)$AUC$cutpoints
    N <- nrow(d)
    pos <- d$score >= 2
    binom.se <- function(p,n){sqrt(p*(1-p)/n*N/(N-1))}
    TPR <- mean(pos[d$Y == 1])
    FPR <- mean(pos[d$Y == 0])
    PPV <- mean(d$Y[pos])
    NPV <- mean(1-d$Y[!pos])
    expect_equal(c(x$TPR,x$FPR,x$PPV,x$NPV),c(TPR,FPR,PPV,NPV))
    expect_equal(c(x$se.TPR,x$se.FPR,x$se.PPV,x$se.NPV),
                 c(binom.se(TPR,sum(d$Y == 1)),binom.se(FPR,sum(d$Y == 0)),
                   binom.se(PPV,sum(pos)),binom.se(NPV,sum(!pos))))
})
# }}}
