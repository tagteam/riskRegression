context("Cause-specific Cox regression")

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





