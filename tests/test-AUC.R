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
    expect_equal(A,B)
    expect_equal(a,b)
})
# }}}

