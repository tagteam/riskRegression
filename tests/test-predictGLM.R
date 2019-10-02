### test-predictGLM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  2 2019 (11:35) 
## Version: 
## Last-Updated: okt  2 2019 (12:13) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Settings
library(riskRegression)
library(testthat)
library(data.table)
library(lava)
context("function predictGLM")


## * [predictGLM] vs. lava
cat("[predictGLM] vs. lava \n")

## ** data
n <- 100
set.seed(10)
dt <- sampleData(n, outcome="binary")

fit <- glm(formula = Y ~ X1+X2, data=dt, family = "binomial")

## ** iid
test_that("[predictGLM] compare to lava",{

    e.RR <- predictGLM(fit, newdata = dt, average.iid = FALSE)
    e.RR.iid <- attr(e.RR,"iid")
    attr(e.RR,"iid") <- NULL

    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        return(expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1")))
    })
    ## check point estimate
    expect_equal(e.RR[,1], predict(fit, newdata = dt, type = "response"))
    expect_equal(unname(e.RR[,1]), unname(e.lava$coef))

    ## check variance
    expect_equal(unname(tcrossprod(e.RR.iid)), unname(e.lava$vcov))
})

## ** average.iid 
test_that("[predictGLM] compare to lava (average.iid, no factor)",{

    ## no factor
    e.RR0 <- predictGLM(fit, newdata = dt, average.iid = FALSE)
    
    e.RR <- predictGLM(fit, newdata = dt, average.iid = TRUE)
    e.RR.iid <- attr(e.RR,"iid")
    attr(e.RR,"iid") <- NULL

    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        eXb <- expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1"))
        return(list(risk = eXb))},
        average=TRUE)

    ## check point estimate
    expect_equal(unname(mean(e.RR[,1])), unname(e.lava$coef))

    ## check iid
    expect_equal(colMeans(attr(e.RR0,"iid")), unname(e.RR.iid[,1]))

    ## check variance
    expect_equal(unname(sum((e.RR.iid + (e.RR-mean(e.RR))[,1]/NROW(e.RR))^2)), unname(e.lava$vcov)[1,1])
})

## ** average.iid with factor
test_that("[predictGLM] compare to lava (average.iid, factor)",{

    factor <- TRUE
    attr(factor,"factor") <- matrix(1:NROW(dt), ncol = 1)
    
    e.RR0 <- predictGLM(fit, newdata = dt, average.iid = FALSE)
    
    e.RR <- predictGLM(fit, newdata = dt, average.iid = factor)
    e.RR.iid <- attr(e.RR,"iid")
    attr(e.RR,"iid") <- NULL
    e.RR[,1] <- e.RR[,1] * (1:NROW(dt))
    
    e.lava <- estimate(fit, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        eXb <- expit(a + b * (data[["X1"]]=="1") + c * (data[["X2"]]=="1"))
        return(list(risk = eXb*(1:NROW(data))))},
        average=TRUE)

    ## check point estimate
    expect_equal(unname(mean(e.RR[,1])), unname(e.lava$coef))

    ## check iid
    expect_equal(colMeans(colMultiply_cpp(attr(e.RR0,"iid"), scale = 1:NROW(dt))), unname(e.RR.iid[,1]))

    ## check variance
    expect_equal(unname(sum((e.RR.iid + (e.RR-mean(e.RR))[,1]/NROW(e.RR))^2)), unname(e.lava$vcov)[1,1])
})

######################################################################
### test-predictGLM.R ends here
