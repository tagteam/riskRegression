context("Prediction error")

test_that("predictions",{
    library(riskRegression)
    library(prodlim)
    data(Melanoma)
    p1 <- predict(ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1),newdata=Melanoma,times=c(0,1,100,1000))
    f2 <- ARR(Hist(time,status)~timevar(sex),data=Melanoma,cause=1)
    p2 <- predict(f2,newdata=Melanoma,times=c(0,1,100,1000))
}

test_that("Brier score",{
    library(riskRegression)
    data(Melanoma)
    fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
    fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
    fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
    library(pec)
    x <- pec(list(ARR=fit.arr2a,ARR.power=fit.arr2b,LRR=fit.lrr),data=Melanoma,formula=Hist(time,status)~1,cause=1,B=10,splitMethod="none")
})
