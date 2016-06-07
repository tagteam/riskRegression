context("Cause-specific Cox regression")


test_that("predictSurv",{
    set.seed(17)
    library(pec)
    library(prodlim)
    d <- prodlim::SimSurv(100)
    f <- coxph(Surv(time,status)~X1+X2,data=d)
    h <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE)
    af <- predictSurvProb(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
    bf <- predictSurv(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
    expect_equal(af,bf)
})

test_that("Cox models",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- CSC(Hist(time,event)~X1+X2,data=d)
    A <- CSC(Hist(time,event)~X1+X2,data=d,survtype="surv")
    a1 <- coxph(Surv(time,event==1)~X1+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~X1+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~X1+X2,data=d)
    expect_equal(coef(a$models[[1]]),coef(a1))
    expect_equal(coef(a$models[[2]]),coef(a2))
    expect_equal(coef(A$models[[2]]),coef(A2))
})

test_that("strata",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- CSC(Hist(time,event)~strata(X1)+X2,data=d)
    A <- CSC(Hist(time,event)~strata(X1)+X2,data=d,survtype="surv")
    a1 <- coxph(Surv(time,event==1)~strata(X1)+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~strata(X1)+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~strata(X1)+X2,data=d)
    expect_equal(coef(a$models[[1]]),coef(a1))
    expect_equal(coef(a$models[[2]]),coef(a2))
    expect_equal(coef(A$models[[2]]),coef(A2))
})

test_that("CSC many character valued causes",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    d$event <- as.character(factor(d$event,labels=c("a","b","c")))
    m1 <- CSC(Hist(time,event)~strata(X1)+X2,data=d)
    m2 <- CSC(Hist(time,event)~strata(X1)+X2,data=d,survtype="surv",cause="b")
    expect_equal(round(coef(m1$models[[2]])[[1]],6),0.535059)
})





