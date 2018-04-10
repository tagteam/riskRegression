library(riskRegression)
library(testthat)
library(survival)
require(rms)
require(survival)
context("Ate checks")
# {{{ header with data  
set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
handler <- if (Sys.info()["sysname"] == "Windows") "foreach" else "mclapply"
verbose=FALSE
# }}}

# {{{ G formula: coxph, cph, sequential, one and several time points
test_that("G formula: coxph, cph, bootstrap sequential, one and several time points",{
    fit.cph <- cph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ate.1a <- ate(fit.cph,data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1,handler=handler, B = 2, y = TRUE, mc.cores=1,verbose=verbose)
    ate.1b <- ate(fit.cph, data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    fit.coxph <- coxph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ate.2a <- ate(fit.coxph,data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    ate.2b <- ate(fit.coxph, data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    attr(ate.2a,"class") <- "NULL"
    attr(ate.2b,"class") <- "NULL"
    attr(ate.1a,"class") <- "NULL"
    attr(ate.1b,"class") <- "NULL"
    expect_equal(ate.1a,ate.2a,tolerance = .00001)
    expect_equal(ate.1b,ate.2b,tolerance = .00001)
})
# }}}
# {{{ Cox model - fully stratified
test_that("G formula: coxph, cph, fully stratified",{
    fit <- cph(formula = Surv(time,event)~ strat(X1),data=dtS,y=TRUE,x=TRUE)
    ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
        times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    # one time point
    fit <- coxph(Surv(time,event)~ strata(X1),data=dtS,y=TRUE,x=TRUE)
    ate(fit,data = dtS, treatment = "X1", contrasts = NULL,
        times = 1, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
})
# }}}
# {{{ Cox model - compare to explicit computation
test_that("Cox model - compare to explicit computation",{
    set.seed(10)
    n <- 5e1
    dtS <- sampleData(n,outcome="survival")
    dtS$time <- round(dtS$time,1)
    dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
    fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ## automatically
    ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                  times = 5:7, B = 0, se = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",lower],
                 c(0.2756147, 0.3220124, 0.3492926),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",upper],
                 c(0.6513235, 0.7011545, 0.7276942),
                 tol = 1e-6)
})
                                        # }}}
                                        # {{{ agreement one timepoint with several timepoints

test_that("Cox model - check internal consistency (one or several timepoints)",{
    ## (previously reported bug)
    d <- sampleData(1000,outcome="survival")
    fit <- coxph(Surv(time,event)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE)
    a1 <- ate(fit,data=d,treatment="X1",time=5,bootci.method="wald")
    a2 <- ate(fit,data=d,treatment="X1",time=5:7,bootci.method="wald")
    expect_equal(a2$riskComparison[time==5,],
                 a1$riskComparison)
    expect_equal(a2$meanRisk[time==5,],
                 a1$meanRisk)
})

# }}}
                                        # {{{ check against manual computation
if (FALSE)
    test_that("check against manual computation",{
        ## automatically
        fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                      times = 5:7, B = 0, se = TRUE, mc.cores=1,handler=handler,
                      verbose=verbose)

        ATE <- list()
        ATE.iid <- list()
        for(iT in c("T0","T1","T2")){
            # iT <- "T0"
            newdata0 <- copy(dtS)
            newdata0$X1 <- iT
            fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
            resPred <- predictCox(fit, newdata = newdata0, time = 5:7, iid = TRUE, log.transform = FALSE)
            ATE[[iT]] <- colMeans(1-resPred$survival)
            ATE.iid_term1 <- apply(-resPred$survival.iid,3,colMeans)
            ATE.iid_term2 <- apply(1-resPred$survival, 1, function(x){x-ATE[[iT]]})/n
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

        expect_equal(diffATE.lower,
                     ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.lower])
        expect_equal(diffATE.upper,
                     ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.upper])

ratioATE <- ATE[["T0"]]/ATE[["T1"]]

ratioATE.iid <- rowScale_cpp(ATE.iid[["T0"]],ATE[["T1"]])-rowMultiply_cpp(ATE.iid[["T1"]], ATE[["T0"]]/ATE[["T1"]]^2)
ratioATE.se <- sqrt(apply(ratioATE.iid^2, 2, sum))
ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se

        expect_equal(ratioATE.lower,
                     ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.lower])
        expect_equal(ratioATE.upper,
                     ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.upper])
})
# }}}
# {{{ CSC model bootstrap
test_that("CSC model bootstrap",{
    df <- sampleData(1e2,outcome="competing.risks")
    df$time <- round(df$time,1)
    df$X1 <- factor(rbinom(1e2, prob = c(0.4,0.3) , size = 2), labels = paste0("T",0:2))
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit, data = df, treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 0, mc.cores=1,handler=handler,verbose=verbose)
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df,  treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 1, mc.cores=1,handler=handler,verbose=verbose)
})
# }}}
# {{{ parallel computation
  
test_that("mcapply vs. foreach",{
    set.seed(10)
    df <- sampleData(3e2,outcome="competing.risks")
    if(parallel::detectCores()>1){
        fit = CSC(formula = Hist(time,event)~ X1+X2, data = df, cause=1)
        time2.mc <- system.time(
            res2.mc <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                           times = 7, cause = 1, B = 2, mc.cores=2, handler = "mclapply", seed = 10, verbose = FALSE)
        )
        time2.for <- system.time(
            res2.for <- ate(fit,data = df, treatment = "X1", contrasts = NULL,
                            times = 7, cause = 1, B = 2, mc.cores=2, handler = "foreach", seed = 10, verbose = FALSE)
        )
        expect_equal(res2.mc, res2.for)
    }})
# }}}
# {{{ LRR
## context("ate for LRR - TO BE DONE")

## n <- 1e2
## dtS <- sampleData(n,outcome="survival")
## dtS$time <- round(dtS$time,1)
## dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

#### NOT GOOD: a different LRR should be fitted at each timepoint
## fitL <- LRR(formula = Hist(time,event)~ X1+X2, data=dtS, cause = 1)
## ate(fitL,data = dtS, treatment = "X1", contrasts = NULL,
##      times = 5:7, B = 3,  mc.cores=1,handler=handler,verbose=verbose)
## ate(fitL, data = dtS, treatment = "X1", contrasts = NULL,
##     times = 7, B = 3, mc.cores=1,handler=handler,verbose=verbose)

# }}}
# {{{ random forest 
## context("ate for random forest")
## n <- 1e2
## dtS <- sampleData(n,outcome="survival")
## dtS$time <- round(dtS$time,1)
## dtS$X1 <- rbinom(n, prob = c(0.3,0.4) , size = 2)
## dtS$X1 <- factor(dtS$X1, labels = paste0("T",unique(dtS$X1)))
## if(require(randomForestSRC)){
## fitRF <- rfsrc(formula = Surv(time,event)~ X1+X2, data=dtS)
## ate(fitRF, data = dtS, treatment = "X1", contrasts = NULL,
## times = 5:7, B = 2, mc.cores=1,handler=handler,verbose=verbose)
## }
# }}}



