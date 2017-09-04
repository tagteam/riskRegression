library(riskRegression)
library(testthat)

#### Cox model ####
context("ate for Cox model")
set.seed(10)

n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

if(require(rms)){
    test_that("G formula: cph, sequential",{
        # one time point
        fit <- cph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
            times = 1, B = 2, y = TRUE, mc.cores=1)
        # several time points
        ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
            times = 1:2, B = 2, y = TRUE, mc.cores=1)
    })
    if(parallel::detectCores()>1){
        test_that("G formula: cph, parallel",{
            fit=cph(formula = Surv(time,event)~ strat(X1)+X2,data=dtS,y=TRUE,x=TRUE)
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
        fit <- coxph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
            times = 1, B = 2, y = TRUE, mc.cores=1)
        # several time points
        ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
            times = 1:2, B = 2, y = TRUE, mc.cores=1)
    })
    if(parallel::detectCores()>1){
        test_that("G formula: coxph, parallel",{
            fit=coxph(formula = Surv(time,event)~ strata(X1)+X2,data=dtS,y=TRUE,x=TRUE)
            set.seed(10)
            ate(object = fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "mclapply", verbose = FALSE)
            set.seed(10)
            ate(object=fit, data = dtS, treatment = "X1", contrasts = NULL,
                times = 1:2, B = 2, y = TRUE, mc.cores=1, handler = "foreach", verbose = FALSE)
        })
    }
}

#### Fully stratified Cox model ####
context("ate for fully stratified Cox model")

n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

if(require(rms)){
  test_that("G formula: cph, fully stratified",{
    fit <- cph(formula = Surv(time,event)~ strat(X1),data=dtS,y=TRUE,x=TRUE)
    ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
        times = 1:2, B = 2, y = TRUE, mc.cores=1)
  })
}

if(require(survival)){
  test_that("G formula: coxph, fully stratified",{
    # one time point
    fit <- coxph(Surv(time,event)~ strata(X1),data=dtS,y=TRUE,x=TRUE)
    ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
        times = 1, B = 2, y = TRUE, mc.cores=1)
    
    predictCox(fit, newdata = dtS, times = 1)
    predictRisk(fit, newdata = dtS, times = 1)
    
  })
}

#### standard error

set.seed(10)
n <- 5e1

dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)

## automatically
ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
              times = 5:7, B = 0, se = TRUE, mc.cores=1)

ateFit

expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",lower],
             c(0.2005312, 0.2467174, 0.2739946),
             tol = 1e-6)
expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",upper],
             c(0.7264070, 0.7764495, 0.8029922),
             tol = 1e-6)

## manually
ATE <- list()
ATE.iid <- list()

for(iT in c("T0","T1","T2")){ # iT <- "T0"
  newdata0 <- copy(dtS)
  newdata0$X1 <- iT
  resPred <- predictCox(fit, newdata = newdata0, time = 5:7, iid = TRUE, logTransform = FALSE)
  ATE[[iT]] <- colMeans(1-resPred$survival)
  
  ATE.iid_term1 <- apply(-resPred$survival.iid,3,colMeans)
  ATE.iid_term2 <- apply(1-resPred$survival, 1, function(x){x-ATE[[iT]]})/sqrt(n)
  ATE.iid[[iT]] <- t(ATE.iid_term1) + t(ATE.iid_term2)
  
  ATE.se <- sqrt(apply(ATE.iid[[iT]]^2, 2, sum))
  ATE.lower <- ATE[[iT]]+ qnorm(0.025) * ATE.se
  ATE.upper <- ATE[[iT]] + qnorm(0.975) * ATE.se
  
  expect_equal(ATE.lower, ateFit$meanRisk[Treatment == iT,lower])
  expect_equal(ATE.upper, ateFit$meanRisk[Treatment == iT,upper])
}

diffATE <- ATE[["T0"]]-ATE[["T1"]]
diffATE.iid <- ATE.iid[["T0"]]-ATE.iid[["T1"]]
diffATE.se <- sqrt(apply(diffATE.iid^2, 2, sum))
diffATE.lower <- diffATE + qnorm(0.025) * diffATE.se
diffATE.upper <- diffATE + qnorm(0.975) * diffATE.se

expect_equal(diffATE.lower, ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.lower])
expect_equal(diffATE.upper, ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.upper])


#### CSC model ####
context("ate for fully CSC model")

df <- sampleData(1e2,outcome="competing.risks")
df$time <- round(df$time,1)
df$X1 <- factor(rbinom(1e2, prob = c(0.4,0.3) , size = 2), labels = paste0("T",0:2))

test_that("no bootstrap",{
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 0, mc.cores=1)
})

test_that("one bootstrap",{
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df,  treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 1, mc.cores=1, verbose = FALSE)
})

#### parallel computation
context("ate with parallel computation")
set.seed(10)
df <- sampleData(3e2,outcome="competing.risks")

fit <- CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
time1.mc <- system.time(
    res1.mc <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                   times = 7, cause = 1, B = 2, mc.cores=1, handler = "mclapply", seed = 10, verbose = FALSE)
)
time1.for <- system.time(
    res1.for <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                    times = 7, cause = 1, B = 2, mc.cores=1, handler = "foreach", seed = 10, verbose = FALSE)
)

res1.mc$meanRisk
res1.for$meanRisk

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

#### LRR ####
context("ate for LRR - TO BE DONE")

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
context("ate for random forest")

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



