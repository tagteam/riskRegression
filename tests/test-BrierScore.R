library(testthat)
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
# {{{ "Brier score"
if (class(try(riskRegression.test,silent=TRUE))[1]!="try-error"){
    test_that("Brier score",{
        library(riskRegression)
        library(survival)
        data(Melanoma)
        fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
        ## predictRisk(fit.lrr,times=c(1,10,100,1000),newdata=Melanoma)
        fit.arr2 <- ARR(Hist(time,status)~thick+age,data=Melanoma,cause=1)
        fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
        fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
        library(pec)
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
