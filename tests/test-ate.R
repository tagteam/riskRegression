library(riskRegression)
library(rms)
library(survival)
library(testthat)
library(data.table)

set.seed(10)
n <- 5e1
dt <- sampleData(n,outcome="competing.risks")
dt[,status:=1*(event!=0)]
dt$time.ties <- round(dt$time)
dt$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
handler <- if (Sys.info()["sysname"] == "Windows") "foreach" else "mclapply"
verbose <- FALSE

set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

test_that("Cox model - compare to explicit computation",{
    fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                  times = 5:7, B = 0, se = TRUE, mc.cores=1,handler=handler,verbose=verbose)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.lower],
                 c(0.2756147, 0.3220124, 0.3492926),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.upper],
                 c(0.6513235, 0.7011545, 0.7276942),
                 tol = 1e-6)
})

test_that("Cox model - check internal consistency (one or several timepoints)",{
    ## (previously reported bug)
    d <- sampleData(1000,outcome="survival")
    fit <- coxph(Surv(time,event)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE)
    a1 <- ate(fit,data=d,treatment="X1",time=5,bootci.method="wald", verbose = verbose)
    a2 <- ate(fit,data=d,treatment="X1",time=5:7,bootci.method="wald", verbose = verbose)

    expect_equal(a2$riskComparison[time==5,],
                 a1$riskComparison)
    expect_equal(a2$meanRisk[time==5,],
                 a1$meanRisk[time==5,])
    })

test_that("check against manual computation",{
    ## automatically
    fit <- cph(formula = Surv(time,event)~ X1+X2,data= dtS,y=TRUE,x=TRUE)
    ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                  times = 5:7,iid=TRUE, B = 0, se = TRUE, mc.cores=1,handler=handler,
                  verbose=verbose)
    ATE <- list()
    ATE.iid <- list()
    for(iT in c("T0","T1","T2")){ ## iT <- "T0"
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
        expect_equal(ATE.lower, ateFit$meanRisk[Treatment == iT,meanRisk.lower])
        expect_equal(ATE.upper, ateFit$meanRisk[Treatment == iT,meanRisk.upper])
    }
    diffATE <- ATE[["T1"]]-ATE[["T0"]]
    diffATE.iid <- ATE.iid[["T1"]]-ATE.iid[["T0"]]
    diffATE.se <- sqrt(apply(diffATE.iid^2, 2, sum))
    diffATE.lower <- diffATE + qnorm(0.025) * diffATE.se
    diffATE.upper <- diffATE + qnorm(0.975) * diffATE.se
    expect_equal(diffATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.lower])
    expect_equal(diffATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.upper])
    ratioATE <- ATE[["T1"]]/ATE[["T0"]]
    ratioATE.iid <- rowScale_cpp(ATE.iid[["T1"]],ATE[["T0"]])-rowMultiply_cpp(ATE.iid[["T0"]], ATE[["T1"]]/ATE[["T0"]]^2)
    ## Brice: I think that this here is how you do it in calcSeATE.R lines 166-168
    ## ratioATE.iid <- ATE.iid[["T1"]]/ATE[["T0"]] -  ATE.iid[["T0"]] * ATE[["T1"]]/ATE[["T0"]]^2
    a <- ATE.iid[["T1"]][,1]
    b <- ateFit$iid$T1[,"5"]
    all.equal(a,b)
                                        # 
    ratioATE.se <- sqrt(rowSums(ratioATE.iid^2))
    ratioATE.se <- sqrt(colSums(ratioATE.iid^2))
    ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
    ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se
    expect_equal(ratioATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.lower])
    expect_equal(ratioATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.upper])
})

test_that("stratified ATE",{
    set.seed(10)
    d <- sampleData(n=100)
    e.coxph <- coxph(Surv(time, event == 1) ~ X1, data = d,
                     x = TRUE, y = TRUE)
    outATE <- ate(e.coxph, data = d, treatment = NULL, strata = "X1", times = 1, verbose = FALSE)

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

## * pre-computation of iidCox does not affect the results
test_that("precompute iid", {
    set.seed(10)
    d <- sampleData(100)

    GS <- ate(Hist(time,event)~X1+X2,
              treatment = "X1", data = d, times = 1:5, cause = 1,
              verbose = FALSE)
    e.CSC <- CSC(Hist(time,event)~X1+X2, data = d)
    e.CSC$models[[1]]$iid <- iidCox(e.CSC$models[[1]])
    e.CSC$models[[2]]$iid <- iidCox(e.CSC$models[[2]])


    test <- ate(e.CSC,
                treatment = "X1", data = d, times = 1:5, cause = 1,
                verbose = FALSE)


    expect_equal(test$meanRisk, GS$meanRisk)
})
