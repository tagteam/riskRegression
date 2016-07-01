library(testthat)
context("G-formula to estimate the average treatment effect")

test_that("no boostrap",{
  res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
         times = 7, cause = 1, B = 0, mc.cores=1)
})

test_that("one boostrap",{
  res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                times = 7, cause = 1, B = 1, mc.cores=1)
})