library(testthat)
context("Cause-specific Cox regression")
library(riskRegression)
library(rms)
library(survival)
library(prodlim)

test_that("CSC vs prodlim",{
    data(Melanoma)
    nd <- data.frame(sex=factor(levels(Melanoma$sex)))
    A <- prodlim(Hist(time,status)~1,data=Melanoma)
    B <- CSC(Hist(time,status)~1,data=Melanoma)
    pB <- predictRisk(B,times=sort(unique(Melanoma$time)),newdata=nd[1,,drop=FALSE],cause=1)
    pA <- predictRisk(A,times=sort(unique(Melanoma$time)),newdata=nd[1,,drop=FALSE],cause=1)
    expect_equal(ignore_attr=TRUE,c(pB),pA,tolerance=1e3)
    a <- prodlim(Hist(time,status)~sex,data=Melanoma)
    b <- CSC(Hist(time,status)~strat(sex),data=Melanoma,fitter="cph")
    c <- CSC(Hist(time,status)~strat(sex),data=Melanoma,surv.type="survival",fitter="cph")
    pa <- predictRisk(a,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    pb <- predictRisk(b,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    pc <- predictRisk(c,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    expect_equal(ignore_attr=TRUE,c(pb),c(pa),tolerance=0.1)
    expect_equal(ignore_attr=TRUE,c(pb),c(pc),tolerance=0.00000001)
#     pd <- predictEventProb(c,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
#     expect_equal(ignore_attr=TRUE,c(pc),c(pd),tolerance=0.00000001)
    u <- CSC(Hist(time,status)~strat(sex)+age+invasion,data=Melanoma,fitter="cph")
    v <- CSC(Hist(time,status)~strat(sex)+age+invasion,data=Melanoma,surv.type="survival",fitter="cph")
    pu <- predictRisk(u,times=c(0,10,100,1000,2000),newdata=Melanoma[c(17,84),],cause=1)
    pv <- predictRisk(v,times=c(0,10,100,1000,2000),newdata=Melanoma[c(17,84),],cause=1)
    expect_equal(ignore_attr=TRUE,c(pu),c(pv),tolerance=0.1)
    ## plot(a,cause=1,ylim=c(0,.5),confint=FALSE)
    ## lines(sort(unique(Melanoma$time)),pb[1,],lty=3)
    ## lines(sort(unique(Melanoma$time)),pb[2,],lty=3,col=2)
})

if (requireNamespace("pec",quietly=TRUE)){
    library(pec)
    test_that("predictSurv",{
        set.seed(17)
        d <- prodlim::SimSurv(100)
        f <- coxph(Surv(time,status)~X1+X2,data=d,x=TRUE,y=TRUE)
        h <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
        af <- predictRisk(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        bf <- 1-predictSurvProb(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        expect_equal(ignore_attr=TRUE,unname(af),unname(bf),tolerance = 1e-8)
        ah <- predictRisk(h,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        bh <- 1-predictSurvProb(h,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        colnames(bh) <- NULL
        expect_equal(ignore_attr=TRUE,unname(ah),unname(bh),tolerance = 1e-8)
    })
}
test_that("Cox models",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- CSC(Hist(time,event)~X1+X2,data=d)
    A <- CSC(Hist(time,event)~X1+X2,data=d,surv.type="surv")
    a1 <- coxph(Surv(time,event==1)~X1+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~X1+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~X1+X2,data=d)
    expect_equal(ignore_attr=TRUE,coef(a$models[[1]]),coef(a1),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(a$models[[2]]),coef(a2),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(A$models[[2]]),coef(A2),tolerance = 1e-8)
})

test_that("strata",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- CSC(Hist(time,event)~strata(X1)+X2,data=d)
    A <- CSC(Hist(time,event)~strata(X1)+X2,data=d,surv.type="surv")
    a1 <- coxph(Surv(time,event==1)~strata(X1)+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~strata(X1)+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~strata(X1)+X2,data=d)
    expect_equal(ignore_attr=TRUE,coef(a$models[[1]]),coef(a1),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(a$models[[2]]),coef(a2),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(A$models[[2]]),coef(A2),tolerance = 1e-8)
})

test_that("strat and strata",{
    data(Melanoma)
    a <- CSC(Hist(time,status)~strat(sex)+age+invasion+logthick+strat(epicel)+strat(ulcer),data=Melanoma,fitter="cph")
    predictRisk(a,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    b <- CSC(Hist(time,status)~strata(sex)+age+invasion+logthick+strata(epicel)+strata(ulcer),data=Melanoma,fitter="coxph")
    pa <- predictRisk(a,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    pb <- predictRisk(b,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    expect_equal(ignore_attr=TRUE,pa,pb,tolerance=1e-6)
})

test_that("factor interactions with strata are identifiable",{
    data(Melanoma)
    Melanoma.numeric <- Melanoma
    Melanoma.numeric$ulcer <- as.integer(Melanoma.numeric$ulcer == "present")
    Melanoma.numeric$epicel <- as.integer(Melanoma.numeric$epicel == "present")
    csc.formula <- list(
        Hist(time,status) ~ sex:(ulcer*age + logthick + epicel) + strata(sex),
        Hist(time,status) ~ sex:(ulcer*age + logthick + epicel) + strata(sex)
    )

    fit.factor <- CSC(csc.formula, data = Melanoma)
    fit.numeric <- CSC(csc.formula, data = Melanoma.numeric)

    expect_true(all(vapply(coef(fit.factor), function(x){all(is.finite(x))}, logical(1))))
    for(iCause in 1:2){
        p.factor <- predictRisk(fit.factor, newdata = Melanoma[1:20,],
                                times = c(100,1000,3000), cause = iCause,
                                product.limit = FALSE)
        p.numeric <- predictRisk(fit.numeric, newdata = Melanoma.numeric[1:20,],
                                 times = c(100,1000,3000), cause = iCause,
                                 product.limit = FALSE)
        expect_equal(ignore_attr = TRUE, p.factor, p.numeric, tolerance = 1e-8)
    }
})

test_that("large stratified predictions keep all covariates when compressing newdata",{
    data(Melanoma)
    Melanoma$ulcer <- as.integer(Melanoma$ulcer == "present")
    Melanoma$epicel <- as.integer(Melanoma$epicel == "present")
    csc.stratified <- CSC(list(
        Hist(time,status) ~ sex:(ulcer*age + logthick + epicel) + strata(sex),
        Hist(time,status) ~ sex:(ulcer*age + logthick + epicel) + strata(sex)
    ), data = Melanoma)
    newdata <- data.frame(
        age = rep(seq(4,95,length.out = 100), length.out = 5000),
        logthick = rep(seq(-2,2,length.out = 101), length.out = 5000),
        sex = factor(rep("Male", 5000), levels = levels(Melanoma$sex)),
        ulcer = rep(0:1, length.out = 5000),
        epicel = rep(0:1, length.out = 5000)
    )

    p.compressed <- predictRisk(csc.stratified, newdata = newdata,
                                times = 4000, cause = 1,
                                product.limit = FALSE)
    p.full <- predictRisk(csc.stratified, newdata = newdata,
                          times = 4000, cause = 1,
                          product.limit = FALSE,
                          store = list(data = "full"))
    expect_equal(ignore_attr = TRUE, p.compressed, p.full, tolerance = 1e-8)
})

# test_that("CSC many character valued causes",{
#     set.seed(17)
#     d <- prodlim::SimCompRisk(100)
#     d$event <- as.character(factor(d$event,labels=c("a","b","c")))
#     m1 <- CSC(Hist(time,event)~strat(X1)+X2,data=d,fitter="cph")
#     m2 <- CSC(Hist(time,event)~strata(X1)+X2,data=d,surv.type="surv",cause="b")
#     expect_equal(ignore_attr=TRUE,round(coef(m1$models[[2]])[[1]],6),0.535059)
# })




