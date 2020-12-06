library(testthat)
library(data.table)
context("Prediction error")
# {{{ "Brier score censored data order"
test_that("Brier score censored data order",{
    library(riskRegression)
    library(survival)
    library(prodlim)
    data(Melanoma)
    setDT(Melanoma)
    fit <- coxph(Surv(time,status!=0)~invasion+epicel+logthick,data=Melanoma,x=TRUE)
    ## many ties in Melanoma
    setkey(Melanoma,age)
    a <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="marginal",metric="Brier")
    A <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="cox",metric="Brier")
    setkey(Melanoma,logthick)
    b <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="marginal",metric="Brier")
    B <- Score(list(fit),data=Melanoma,Surv(time,status!=0)~invasion+epicel+logthick,cens.model="cox",metric="Brier")
    a$call <- b$call <- A$call <- B$call <- NULL
    ## expect_error(expect_equal(a,b,tolerance = .002))
    expect_equal(a,b,tolerance = .02)
    expect_equal(A,B,tolerance=.02)
})

# }}}
