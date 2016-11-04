library(riskRegression)
library(testthat)
context("G-formula to estimate the average treatment effect")

#### Cox model ####
n <- 1e2
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

if(require(rms)){
    test_that("G formula: cph, sequential",{
        # one time point
        fit <- cph(Surv(time,event)~ X1+X2,data=dtS,x=TRUE,surv=TRUE,y=TRUE)
        ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
            times = 1, B = 2, y = TRUE, mc.cores=1)
        # several time points
        ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
            times = 1:10, B = 2, y = TRUE, mc.cores=1)
    })
    if(parallel::detectCores()>1){
        test_that("G formula: cph, parallel",{
            fit=cph(formula = Surv(time,event)~ strat(X1)+X2,data=dtS,surv=TRUE,y=TRUE)
            ate(object = fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "mclapply", verbose = FALSE)
            ate(object=fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "foreach", verbose = FALSE)
        })
    }
}

if(require(survival)){
    test_that("G formula: coxph, sequential",{
        # one time point
        fit <- coxph(Surv(time,event)~ X1+X2,data=dtS,x=TRUE)
        ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
            times = 1, B = 2, y = TRUE, mc.cores=1)
        # several time points
        ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
            times = 1:10, B = 2, y = TRUE, mc.cores=1)
    })
    if(parallel::detectCores()>1){
        test_that("G formula: coxph, parallel",{
            fit=coxph(formula = Surv(time,event)~ strata(X1)+X2,data=dtS)
            ate(object = fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "mclapply", verbose = FALSE)
            ate(object=fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "foreach", verbose = FALSE)
        })
    }
}

#### Fully stratified Cox model ####
n <- 1e2
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

if(require(rms)){
  test_that("G formula: cph, fully stratified",{
    fit <- cph(formula = Surv(time,event)~ strat(X1),data=dtS,y=TRUE)
    ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
        times = 1:10, B = 2, y = TRUE, mc.cores=1)
  })
}

if(require(survival)){
  test_that("G formula: coxph, fully stratified",{
    # one time point
    fit <- coxph(Surv(time,event)~ strata(X1),data=dtS,x=TRUE)
    ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
        times = 1, B = 2, y = TRUE, mc.cores=1)
    
    predictCox(fit, newdata = dtS, times = 1)
    predictRisk(fit, newdata = dtS, times = 1)
    
  })
}

#### CSC model ####
df <- sampleData(1e2,outcome="competing.risks")
df$time <- round(df$time,1)
df$X1 <- factor(rbinom(1e2, prob = c(0.4,0.3) , size = 2), labels = paste0("T",0:2))

test_that("no boostrap",{
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 0, mc.cores=1)
})

test_that("one boostrap",{
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df,  treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 1, mc.cores=1, verbose = FALSE)
})

#### parallel computation
fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
time1.mc <- system.time(
    res1.mc <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                   times = 7, cause = 1, B = 2, mc.cores=1, handler = "mclapply", seed = 10, verbose = FALSE)
)
time1.for <- system.time(
    res1.for <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                    times = 7, cause = 1, B = 2, mc.cores=1, handler = "foreach", seed = 10, verbose = FALSE)
)

test_that("mcapply vs. foreach",{
  expect_equal(res1.mc,res1.for)
})


if(parallel::detectCores()>1){
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    time2.mc <- system.time(
        res2.mc <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                       times = 7, cause = 1, B = 2, mc.cores=2, handler = "mclapply", seed = 10, verbose = FALSE)
    )
    time2.for <- system.time(
        res2.for <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                        times = 7, cause = 1, B = 2, mc.cores=2, handler = "foreach", seed = 10, verbose = FALSE)
    )
  
  test_that("mcapply vs. foreach",{
    expect_equal(res2.mc, res2.for)
  })
  
}


#### Strata with different event times
#boxplot(time ~ X2, df)
df <- sampleData(5e2,outcome="competing.risks")
df$timeC <- df$time
df$eventC <- df$event
df[(df$X2 == 1)*(df$time>7)>0, "timeC"] <- 7
df[(df$X2 == 1)*(df$time>7)>0, "timeC"] <- 0

# plot(prodlim(Hist(timeC,eventC)~ X2+X1, df))

test_that("early strata stop",{
  
  # df[,max(timeC),by = c("X1","X2")]
  
  fit <- CSC(formula =  Hist(timeC,eventC)~ strat(X1) + strat(X2), data = df,cause=1)
  
  res <- ate(fit, data = df, treatment = "X1", contrasts = NULL,
             times = 6.8, cause = 1, B = 1, mc.cores=1, seed = 1, verbose = FALSE)
})

#### LRR ####
n <- 1e2
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))


#### NOT GOOD: a different LRR should be fitted at each timepoint
# fitL <- LRR(formula = Hist(time,event)~ X1+X2, data=dtS)
# ate(fitL,data = dtS, treatment = "X1", contrasts = NULL,
#     times = 5:7, B = 3,  mc.cores=1)
# ate(fitL, data = dtS, treatment = "X1", contrasts = NULL,
#     times = 7, B = 3, mc.cores=1)

#### random forest ####
n <- 1e2
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- rbinom(n, prob = c(0.3,0.4) , size = 2)
dtS$X1 <- factor(dtS$X1, labels = paste0("T",unique(dtS$X1)))

if(require(randomForestSRC)){
    fitRF <- rfsrc(formula = Surv(time,event)~ X1+X2, data=dtS)
    ate(fitRF, data = dtS, treatment = "X1", contrasts = NULL,
        times = 5:7, B = 2, mc.cores=1)
}

