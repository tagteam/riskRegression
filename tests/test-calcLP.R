library(riskRegression)
library(testthat)
library(rms)
library(survival)
context("Compute the linear predictor of a Cox model")

d <- sampleData(1e2, outcome = "survival")

##
test_that("linear predictor (one variable)",{
  mCox <- cph(Surv(time, event) ~ X1, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1"), stratavars = NULL)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(lpTerms, lpGS)
  
  mCox <- coxph(Surv(time, event) ~ X1, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1"), stratavars = NULL)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(unname(lpTerms), lpGS)
})

test_that("linear predictor (two variables)",{
  mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(unname(lpTerms), lpGS)
  
  mCox <- cph(Surv(time, event) ~ X1+X2, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(lpTerms, lpGS)
})

test_that("linear predictor (several variables + interactions)",{
  mCox <- coxph(Surv(time, event) ~ X1*X2+X5+X6*X3, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1","X2","X3","X5","X6"), stratavars = NULL)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(unname(lpTerms), lpGS)
  
  mCox <- cph(Surv(time, event) ~ X1*X2+X5+X6*X3, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = c("X1","X2","X3","X5","X6"), stratavars = NULL, aggregate = TRUE)
  
  lpGS <- predict(mCox, data = d, type = "lp")
  expect_equal(lpTerms, lpGS)
})
             

####
if(FALSE){ # problem with rms when using continous variables
  
  dsave <- d
  d$X6 <- scale(d$X6)
  
  variable <- "X6"
  f <- eval(parse(text = paste0(
    "Surv(time, event) ~ ",variable
  )))
  
  mCox1 <- coxph(f, data = d)
  lpTerms1 <- lpCox(mCox1, data = d, xvars = variable, stratavars = NULL, aggregate = FALSE, centered = TRUE)
  
  mCox <- cph(f, data = d)
  lpTerms <- lpCox(mCox, data = d, xvars = variable, stratavars = NULL, aggregate = FALSE, centered = TRUE)
  
  range(predict(mCox1, d, type = "lp") - lpTerms1)
  range(predict(mCox, d, type = "lp") - lpTerms)
  
  
  range(predict(mCox, d, type = "lp") - predict(mCox1, d, type = "lp"))
  range(predict(mCox, d, type = "terms") - predict(mCox1, d, type = "terms"))
  
  
  Xb_terms <- rms:::predictrms(mCox, d, type = "terms", se.fit = FALSE, conf.int = FALSE, conf.type ="mean", 
                               kint = 1, na.action = na.keep, expand.na = TRUE, center.terms = TRUE)
  Xb_lp <- rms:::predictrms(mCox, d, type = "lp", se.fit = FALSE, conf.int = FALSE, conf.type ="mean", 
                            kint = 1, na.action = na.keep, expand.na = TRUE, center.terms = FALSE)
  
  
}






##





