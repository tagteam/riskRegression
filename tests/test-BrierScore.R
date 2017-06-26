library(testthat)
context("Prediction error")

test_that("predictions",{
    library(riskRegression)
    library(survival)
    library(prodlim)
    data(Melanoma)
    p1 <- stats::predict(ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1),newdata=Melanoma,times=c(0,1,100,1000))
    f2 <- ARR(Hist(time,status)~timevar(sex),data=Melanoma,cause=1)
    p2 <- stats::predict(f2,newdata=Melanoma,times=c(0,1,100,1000))
})

test_that("Brier score",{
    library(riskRegression)
    library(survival)
    data(Melanoma)
    fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
    ## predictRisk(fit.lrr,times=c(1,10,100,1000),newdata=Melanoma)
    fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
    fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
    library(pec)
     old <- pec(list(ARR=fit.arr2a,ARR.power=fit.arr2b,LRR=fit.lrr),
               data=Melanoma,
               formula=Hist(time,status)~1,
               cause=1, B=10,splitMethod="none")
    new <- Score(list(ARR=fit.arr2a,ARR.power=fit.arr2b,LRR=fit.lrr),
                 data=Melanoma,
                 times=c(0,sort(unique(Melanoma$time))),
                 metrics="brier",plots=NULL,summary=NULL,
                 formula=Hist(time,status)~1,
                 cause=1, B=10,splitMethod="none")
    nix <- lapply(1:4,function(m){
        expect_equal(new$Brier$score[model==names(new$models)[m]][["Brier"]],old$AppErr[[names(old$AppErr)[[m]]]])})

})
