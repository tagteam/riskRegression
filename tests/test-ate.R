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
verbose <- FALSE
                                        # }}}
                                        # {{{ G formula: coxph, cph, sequential, one and several time points

test_that("G formula: coxph, cph, bootstrap sequential, one and several time points",{
    fit.cph <- cph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ate.1a <- ate(fit.cph,data = dtS, treatment = "X1", contrasts = NULL,seed=3, 
                  times = 1, handler=handler, B = 2, y = TRUE, mc.cores=1,verbose=verbose)
    ate.1a <- confint(ate.1a, bootci.method = "quantile")
    ate.1b <- ate(fit.cph, data = dtS, treatment = "X1", contrasts = NULL,seed=3, 
                  times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    ate.1b <- confint(ate.1b, bootci.method = "quantile")

    fit.coxph <- coxph(Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ate.2a <- ate(fit.coxph,data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    ate.2a <- confint(ate.2a, bootci.method = "quantile")
    ate.2b <- ate(fit.coxph, data = dtS, treatment = "X1", contrasts = NULL,seed=3,
                  times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    ate.2b <- confint(ate.2b, bootci.method = "quantile")
    
    attr(ate.2a,"class") <- "NULL"
    attr(ate.2b,"class") <- "NULL"
    attr(ate.1a,"class") <- "NULL"
    attr(ate.1b,"class") <- "NULL"
    expect_equal(ate.1a,ate.2a,tolerance = .00001)
    expect_equal(ate.1b,ate.2b,tolerance = .00001)

    ## quantile(ate.1a$boot$t[,1], 0.025) 
    expect_equal(ate.1a$meanRisk$meanRiskBoot,
                 c(0.04778180, 0.03383714, 0.03132013), tol = 1e-6)
    expect_equal(ate.1a$meanRisk$se,
                 c(0.02308268, 0.03008728, 0.02586185), tol = 1e-6)
    expect_equal(ate.1a$meanRisk$lower,
                 c(0.03227598, 0.01362597, 0.01394739), tol = 1e-6)
    expect_equal(ate.1a$meanRisk$upper,
                 c(0.06328763, 0.05404831, 0.04869286), tol = 1e-6)
    
    expect_equal(ate.1a$riskComparison$diffMeanBoot,
                 c(0.013944664, 0.016461677, 0.002517013), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$diff.se,
                 c(0.007004595, 0.002779170, 0.004225425), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$diff.lower,
                 c(0.0092393173, 0.0145947658, -0.0003214227), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$diff.upper,
                 c(0.018650011, 0.018328589, 0.005355449), tol = 1e-6)
    ## expect_equal(ate.1a$riskComparison$diff.p.value,
    ## c(4.650420e-02, 3.156692e-09, 5.513872e-01), tol = 1e-6)
  
    expect_equal(ate.1a$riskComparison$ratioMeanBoot,
                 c(1.833739, 1.853041, 1.037422), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$ratio.se,
                 c(0.9483519, 0.7931128, 0.1040106), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$ratio.lower,
                 c(1.1966818, 1.3202665, 0.9675527), tol = 1e-6)
    expect_equal(ate.1a$riskComparison$ratio.upper,
                 c(2.470795, 2.385816, 1.107291), tol = 1e-6)
    ## expect_equal(ate.1a$riskComparison$ratio.p.value,
                 ## c(5.316164e-02, 1.946959e-02, 1.976976e-23), tol = 1e-6)

})

                                        # }}}
                                        # {{{ Cox model - fully stratified

test_that("G formula: coxph, cph, fully stratified",{
    fit <- cph(formula = Surv(time,event)~ strat(X1),data=dtS,y=TRUE,x=TRUE)
    ate2 <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL, seed = 3,
                times = 1:2, B = 2, y = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    ## one time point
    fit <- coxph(Surv(time,event)~ strata(X1),data=dtS,y=TRUE,x=TRUE)
    ate1 <- ate(fit,data = dtS, treatment = "X1", contrasts = NULL, seed = 3,
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
    ateFit <- confint(ateFit)
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
                 a1$meanRisk[time==5,])
})

                                        # }}}
                                        # {{{ check against manual computation
test_that("check against manual computation",{
    ## automatically
    fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                  times = 5:7, B = 0, se = TRUE, mc.cores=1,handler=handler,
                  verbose=verbose)
    ateFit <- confint(ateFit)
    
    ATE <- list()
    ATE.iid <- list()
    for(iT in c("T0","T1","T2")){
                                        # iT <- "T0"
        newdata0 <- copy(dtS)
        newdata0$X1 <- iT
        fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        resPred <- predictCox(fit, newdata = newdata0, time = 5:7, iid = TRUE)
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
    Sres <- print(res, trace = FALSE)
    fit=CSC(formula = Hist(time,event)~ X1+X2, data = df,cause=1)
    res <- ate(fit,data = df,  treatment = "X1", contrasts = NULL,
               times = 7, cause = 1, B = 2, mc.cores=1,handler=handler,verbose=verbose)
    Sres <- print(res, trace = FALSE)
})
                                        # }}}

                                        # {{{ average risk

test_that("stratified ATE",{
    set.seed(10)
    d <- sampleData(n=100)

    e.coxph <- coxph(Surv(time, event == 1) ~ X1, data = d,
                     x = TRUE, y = TRUE)

    outATE <- ate(e.coxph, data = d, treatment = NULL, strata = "X1", time = 1)

    outPred <- predictCox(e.coxph,
                          newdata = d,
                          times = 1,
                          se = TRUE,
                          keep.newdata = TRUE,
                          type = "survival")
    GS.se <- as.data.table(outPred)[,.SD[1],by = "X1"]
    setkeyv(GS.se, cols = "X1")
    test.se <- outATE$meanRisk
    expect_equal(1-GS.se$survival,test.se$meanRisk)
    expect_equal(GS.se$survival.se,GS.se$survival.se)
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



