library(riskRegression)
library(testthat)
context("G-formula to estimate the average treatment effect")

dt <- sampleData(1e3,outcome="competing.risks")
dt$time <- round(dt$time,1)
seqtimes <- sample(x = unique(sort(dt$time)), size = 100) 
dt$X1 <- factor(rbinom(1e3, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))

test_that("no boostrap",{
  res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
         times = 7, cause = 1, B = 0, mc.cores=1)
})

test_that("one boostrap",{
  res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                times = 7, cause = 1, B = 1, mc.cores=1)
})

#### parallel computation

time1.mc <- system.time(
  res1.mc <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                times = 7, cause = 1, B = 50, mc.cores=1, handler = "mclapply", seed = 10)
)
time1.for <- system.time(
  res1.for <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                times = 7, cause = 1, B = 50, mc.cores=1, handler = "foreach", seed = 10)
)

test_that("mcapply vs. foreach",{
  expect_equal(res1.mc,res1.for)
})


if(parallel::detectCores()>1){
  time2.mc <- system.time(
    res2.mc <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                  times = 7, cause = 1, B = 50, mc.cores=2, handler = "mclapply", seed = 10)
  )
  time2.for <- system.time(
    res2.for <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
                  times = 7, cause = 1, B = 50, mc.cores=2, handler = "foreach", seed = 10)
  )
  
  test_that("mcapply vs. foreach",{
    expect_equal(res2.mc, res2.for)
  })
  
}
# system.time(
#   res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
#                 times = 7, cause = 1, B = 10, mc.cores=1, handler = "foreach")
# )
# system.time(
#   res <- ateCSC(formula = Hist(time,event)~ X1+X2,data = dt, treatment = "X1", contrasts = NULL,
#                 times = 7, cause = 1, B = 10, mc.cores=4, handler = "foreach")
# )
