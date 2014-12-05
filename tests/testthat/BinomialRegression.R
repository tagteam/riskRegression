context("Binomial regression")

test_that("Absolute risk regression",{
    set.seed(17)
    library(riskRegression)
    library(prodlim)
    d <- prodlim::SimCompRisk(100)
    a <- ARR(Hist(time,event)~X1+X2,data=d)
    const <- function(x)x
    b <- timereg::comp.risk(Hist(time,event)~ const(X1)+const(X2),data=d,cause=1,model="rcif")
    expect_equal(as.numeric(a$timeConstantEffects$coef),c(b$gamma))
    A <- ARR(Hist(time,event)~X1+strata(X2),data=d)
    B <- timereg::comp.risk(Hist(time,event)~ const(X1)+X2,data=d,cause=1,model="rcif")
    B <- timereg::comp.risk(Event(time,event)~ const(X1)+X2,data=d,cause=1,model="rcif")
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})
test_that("Logistic risk regression",{
    set.seed(17)
    d <- prodlim::SimCompRisk(100)
    a <- LRR(Hist(time,event)~X1+X2,data=d)
    const <- function(x)x
    b <- timereg::comp.risk(Hist(time,event)~ const(X1)+const(X2),data=d,cause=1,model="logistic")
    expect_equal(as.numeric(a$timeConstantEffects$coef),c(b$gamma))
    A <- LRR(Hist(time,event)~X1+strata(X2),data=d)
    B <- timereg::comp.risk(Hist(time,event)~ const(X1)+X2,data=d,cause=1,model="logistic")
    ## head(A$timeVaryingEffects$coef)
    ## head(B$cum)
    expect_equal(as.numeric(A$timeConstantEffects$coef),c(B$gamma))    
})

test_that("Censoring model",{
    library(riskRegression)
    data(Melanoma)
    f1 <- ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox")
    f1a <- ARR(Hist(time,status)~thick+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~thick+strat(invasion)+epicel)
    f2 <- ARR(Hist(time,status)~sex+strata(invasion)+epicel,data=Melanoma,cens.model="cox",cens.formula=~logthick+age)
})


