library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(testthat)
library(data.table)
library(ipw)

## * [ate] G-formula only
cat("[ate] G-formula only \n")

## ** data
set.seed(10)
n <- 5e1
dt <- sampleData(n,outcome="competing.risks")
dt[,status:=1*(event!=0)]
dt$time.ties <- round(dt$time)
dt$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
handler <- if (Sys.info()["sysname"] == "Windows") "foreach" else "mclapply"
verbose <- FALSE

## stratified version
set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

## ** Compare to fix value / explicit computation 
test_that("[ate] Cox model - compare to fix values / explicit computation",{
    seqTime <- c(1,5,10)
    
    fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
    
    ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
                  times = seqTime, B = 0, iid = TRUE, se = TRUE, verbose = verbose)

    ## point estimate
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk],
                 c(0.05842145, 0.46346908, 0.66107744),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.lower],
                 c(0.0000000, 0.2756147, 0.4670461),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.upper],
                 c(0.1178148, 0.6513235, 0.8551088),
                 tol = 1e-6)

    ## standard error
    ATE <- list()
    ATE.iid <- list()
    for(iT in c("T0","T1","T2")){ ## iT <- "T0"
        newdata0 <- copy(dtS)
        newdata0$X1 <- iT
        fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        resPred <- predictCox(fit, newdata = newdata0, time = seqTime, iid = TRUE)
        ATE[[iT]] <- colMeans(1-resPred$survival)
        ATE.iid_term1 <- apply(-resPred$survival.iid,3,colMeans)
        ATE.iid_term2 <- apply(1-resPred$survival, 1, function(x){x-ATE[[iT]]})/n
        ATE.iid[[iT]] <- t(ATE.iid_term1) + t(ATE.iid_term2)
        ATE.se <- sqrt(apply(ATE.iid[[iT]]^2, 2, sum))
        ATE.lower <- ATE[[iT]]+ qnorm(0.025) * ATE.se
        ATE.upper <- ATE[[iT]] + qnorm(0.975) * ATE.se
        expect_equal(pmax(0,ATE.lower), ateFit$meanRisk[Treatment == iT,meanRisk.lower])
        expect_equal(pmin(1,ATE.upper), ateFit$meanRisk[Treatment == iT,meanRisk.upper])
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
    expect_equal(unname(ATE.iid[["T1"]]),unname(ateFit$iid$T1))

    ratioATE.se <- sqrt(rowSums(ratioATE.iid^2))
    ratioATE.se <- sqrt(colSums(ratioATE.iid^2))
    ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
    ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se
    expect_equal(ratioATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.lower])
    expect_equal(ratioATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.upper])
})

test_that("[ate] Cox model - compare to fix values / explicit computation",{
    seqTime <- c(1,5,10)
    
    fitS <- cph(formula = Surv(time,event)~ strat(X1)+X2*X6,data=dtS,y=TRUE,x=TRUE)
    
    ateFit <- ate(fitS, data = dtS, treatment = "X1", contrasts = NULL,
                  times = seqTime, B = 0, iid = TRUE, se = TRUE, verbose = verbose)

    ## point estimate
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk],
                 c(0.0000000, 0.4124494, 0.6420487),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.lower],
                 c(0.0000000, 0.2446007, 0.5004145),
                 tol = 1e-6)
    expect_equal(ateFit$meanRisk[ateFit$meanRisk$Treatment == "T0",meanRisk.upper],
                 c(0.0000000, 0.5802981, 0.7836830),
                 tol = 1e-6)

    ## standard error
    ATE <- list()
    ATE.iid <- list()
    for(iT in c("T0","T1","T2")){ ## iT <- "T0"
        newdata0 <- copy(dtS)
        newdata0$X1 <- iT
        fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
        resPred <- predictCox(fitS, newdata = newdata0, time = seqTime, iid = TRUE)
        ATE[[iT]] <- colMeans(1-resPred$survival)
        ATE.iid_term1 <- apply(-resPred$survival.iid,3,colMeans)
        ATE.iid_term2 <- apply(1-resPred$survival, 1, function(x){x-ATE[[iT]]})/n
        ATE.iid[[iT]] <- t(ATE.iid_term1) + t(ATE.iid_term2)
        ATE.se <- sqrt(apply(ATE.iid[[iT]]^2, 2, sum))
        ATE.lower <- ATE[[iT]]+ qnorm(0.025) * ATE.se
        ATE.upper <- ATE[[iT]] + qnorm(0.975) * ATE.se
        expect_equal(pmax(0,ATE.lower), ateFit$meanRisk[Treatment == iT,meanRisk.lower])
        expect_equal(pmin(1,ATE.upper), ateFit$meanRisk[Treatment == iT,meanRisk.upper])
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
    expect_equal(unname(ATE.iid[["T1"]]),unname(ateFit$iid$T1))

    ratioATE.se <- sqrt(rowSums(ratioATE.iid^2))
    ratioATE.se <- sqrt(colSums(ratioATE.iid^2))
    ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
    ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se
    expect_equal(ratioATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.lower])
    expect_equal(ratioATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.upper])
})


## * [ate] Survival case
cat("[ate] survival case")
## ** no censoring - manual computation
n <- 5e1
tau <- 1.5
## simulate data
set.seed(10)
dtS <- sampleData(n,outcome="survival")
dtS$event <- 1
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$X1 <- as.numeric(dtS$X1)-1

## fit
test_that("check vs. manual calculations", {
    e.S <- coxph(Surv(time, event) ~ X1 + X2 + X3,data = dtS,x = TRUE)
    e.T <- glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")) ## dtS$X1
    e.ateRobust <- ateRobust(data = dtS, times = tau,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.treatment = X1 ~ 1, se = TRUE,
                             type = "survival")

    ## point estimate
    dtS1 <- data.table(X1 = 1, X2 = dtS$X2, X3 = dtS$X3)
    pred.logit <- riskRegression:::predictGLM(e.T, newdata = dtS1, average.iid = FALSE)
    iid.logit <- attr(pred.logit, "iid")
    attr(pred.logit, "iid") <- NULL
    pred.Surv <- predictCox(e.S, newdata = dtS1, times = tau, type = "survival", iid = TRUE)
    iid.Surv <- -pred.Surv$survival.iid[,1,]
    dtS[, pi := as.double(pred.logit)]
    dtS[, r := as.double(1-pred.Surv$survival)]
    expect_equal(dtS[,mean(r)],e.ateRobust$ate.value["risk.1","Gformula"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.IPCW"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.AIPCW"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.AIPCW_knownNuisance"])
  

    ## standard error
    iid.Gformula1 <- (dtS$r - mean(dtS$r))/NROW(dtS)
    iid.Gformula2 <- colMeans(iid.Surv)
    iid.Gformula <- iid.Gformula1 + iid.Gformula2
    expect_equal(sqrt(sum(iid.Gformula^2)),e.ateRobust$ate.se["risk.1","Gformula"])
    iid.IPTW1 <- (dtS[,X1*Y/pi] - dtS[,mean(X1*Y/pi)])/NROW(dtS)
    iid.IPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*dtS$Y/dtS$pi^2) )
    iid.IPTW <- iid.IPTW1 + iid.IPTW2
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.IPCW"])
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.AIPCW"])
    iid.AIPTW1 <- (dtS[,X1*Y/pi + r*(1-X1/pi)] - dtS[,mean(X1*Y/pi + r*(1-X1/pi))])/NROW(dtS)
    iid.AIPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*(dtS$Y-dtS$r)/dtS$pi^2) )
    iid.AIPTW3 <- colMeans( riskRegression::colMultiply_cpp(iid.Surv, scale = (1-dtS$X1/dtS$pi)) )
    iid.AIPTW <- iid.AIPTW1 + iid.AIPTW2 + iid.AIPTW3
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.AIPCW_knownNuisance"])
    ## (sqrt(sum(iid.AIPTW^2)) - sqrt(sum(iid.AIPTW1^2))) / sqrt(sum(iid.AIPTW^2))

    ## compare to ipw package
    e2.ateRobust <- ateRobust(data = dtS, times = tau,
                              formula.event = Surv(time,event) ~ X1 + X2 + X3,
                              formula.treatment = X1 ~ X2+X3, se = TRUE,
                              type = "survival")
    RR.weight <- as.matrix(1/e2.ateRobust$weight.treatment)
    
    ww <- ipwpoint(exposure=X1,family="binomial",link="logit",
                   denominator=~X2+X3,
                   data=dtS)$ipw.weights
    expect_equal(ww, RR.weight[1:NROW(RR.weight) + dtS$X1 * NROW(RR.weight)])

    ## compare with ate
    dtS$X1f <- as.factor(dtS$X1)
    e.ateRR <- ate(event = cph(Surv(time,event) ~ X1f + X2 + X3, data = dtS, x = TRUE, y = TRUE),
                   treatment = glm(X1f ~ 1, data = dtS, family = binomial(link = "logit")),
                   data = dtS, times = tau, verbose = 0,
                   )
    e.ateIP <- ate(event = c("time","event"),
                   treatment = glm(X1f ~ 1, data = dtS, family = binomial(link = "logit")),
                   data = dtS, times = tau, verbose = 0
                   )

    expect_equal(as.double(e.ateRobust$ate.value[,"AIPTW.AIPCW_estimatedNuisance"]),
                 c(e.ateRR$meanRisk[,meanRisk],e.ateRR$riskComparison[,diff]),
                 tol = 1e-4)
    expect_equal(as.double(e.ateRobust$ate.se[,"AIPTW.AIPCW_estimatedNuisance"]),
                 c(e.ateRR$meanRisk[,meanRisk.se],e.ateRR$riskComparison[,diff.se]),
                 tol = 1e-4)

    expect_equal(as.double(e.ateRobust$ate.value[,"IPTW.IPCW"]),
                 c(e.ateIP$meanRisk[,meanRisk],e.ateIP$riskComparison[,diff]),
                 tol = 1e-4)
    expect_equal(as.double(e.ateRobust$ate.se[,"IPTW.IPCW"]),
                 c(e.ateIP$meanRisk[,meanRisk.se],e.ateIP$riskComparison[,diff.se]),
                 tol = 1e-4)
})

## ** censoring - agreement with ate
set.seed(10)
n <- 5e1
tau <- 3
dtS <- sampleData(n,outcome="survival")
e.cox <- coxph(Surv(time, event) ~ X1 + X2 + X3,
               data = dtS,
               x = TRUE)

test_that("Agreement ate-ateRobust (survival)",{
    ## dtS[event == 0, min(time)]
    
    e.ateG <- ate(e.cox, treatment = "X1", times = tau, data = dtS, se = TRUE, verbose = 0)
    e.ateRR <- ate(event = cph(Surv(time,event) ~ X1 + X2 + X3, data = dtS, x = TRUE, y = TRUE),
                   treatment = glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")),
                   censor = cph(Surv(time,event==0) ~ X1, data = dtS, x = TRUE, y = TRUE),
                   data = dtS, times = tau, verbose = 0
                   )
    e.ateIP <- ate(event = c("time","event"),
                   treatment = glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")),
                   censor = cph(Surv(time,event==0) ~ X1, data = dtS, x = TRUE, y = TRUE),
                   data = dtS, times = tau, verbose = 0
                   )

    e.ateRobust <- ateRobust(data = dtS, times = tau,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             se = TRUE, product.limit = FALSE,
                             type = "survival")

    e.ateRobustPL <- ateRobust(data = dtS, times = tau,
                               formula.event = Surv(time,event) ~ X1 + X2 + X3,
                               formula.censor = Surv(time,event==0) ~ X1,
                               formula.treatment = X1 ~ 1,
                               se = TRUE, product.limit = TRUE,
                               type = "survival")
    
    ## ateRobust with product.limit = FALSE agree with ate
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula"]),
                 c(e.ateG$meanRisk[,meanRisk],e.ateG$riskComparison[,diff]),
                 tol = 1e-4)
    expect_equal(as.double(e.ateRobust$ate.se[,"Gformula"]),
                 c(e.ateG$meanRisk[,meanRisk.se],e.ateG$riskComparison[,diff.se]),
                 tol = 1e-4)

    expect_equal(as.double(e.ateRobust$ate.value[,"IPTW.IPCW"]),
                 c(e.ateIP$meanRisk[,meanRisk],e.ateIP$riskComparison[,diff]),
                 tol = 1e-4)
    ## missing term for the influence function in ateRobust
    ## expect_equal(as.double(e.ateRobust$ate.se[,"IPTW.IPCW"]),
                 ## c(e.ateIP$meanRisk[,meanRisk.se],e.ateIP$riskComparison[,diff.se]),
                 ## tol = 1e-4)

    expect_equal(as.double(e.ateRobust$ate.value[,"AIPTW.AIPCW_estimatedNuisance"]),
                 c(e.ateRR$meanRisk[,meanRisk],e.ateRR$riskComparison[,diff]),
                 tol = 1e-4)
    ## missing term for the influence function in ateRobust
    ## expect_equal(as.double(e.ateRobust$ate.se[,"AIPTW.AIPCW_estimatedNuisance"]),
                 ## c(e.ateRR$meanRisk[,meanRisk.se],e.ateRR$riskComparison[,diff.se]),
                 ## tol = 1e-4)

    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.19447, 0.54518, 0.35071), 
                  "IPTW.IPCW" = c(0.15098, 0.71242, 0.56145), 
                  "AIPTW.IPCW_knownNuisance" = c(0.1463, 0.75177, 0.60547), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.1463, 0.75177, 0.60547), 
                  "IPTW.AIPCW" = c(0.15097, 0.70596, 0.55499), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.14629, 0.7453, 0.59901),
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.14629, 0.7453, 0.59901))
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)

    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.04724, 0.13442, 0.12753), 
                  "IPTW.IPCW" = c(0.05683, 0.14767, 0.15823), 
                  "AIPTW.IPCW_knownNuisance" = c(0.05662, 0.15129, 0.15937), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.05571, 0.14053, 0.14949), 
                  "IPTW.AIPCW" = c(0.05681, 0.14873, 0.15921), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.05661, 0.15154, 0.15965), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0557, 0.14153, 0.15041)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
})

## * [ate] Competing risk case
cat("[ate] Competing risk case \n")

## ** no censoring - manual computation
n <- 1e2
tau <- 1.5

## simulate data
set.seed(10)
dtS <- sampleData(n,outcome="competing.risks")
dtS <- dtS[event!=0]
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$X1 <- as.numeric(dtS$X1)-1

## fit
e.S <- CSC(Hist(time, event) ~ X1 + X2 + X3,
           data = dtS)
e.T <- glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")) ## dtS$X1

e.ateRobust <- ateRobust(data = dtS, times = tau,
                         formula.event = Hist(time,event) ~ X1 + X2 + X3,
                         formula.censor = Hist(time,event==0) ~ X1,
                         formula.treatment = X1 ~ 1, se = TRUE, cause = 1,
                         type = "competing.risks")

## predict
dtS1 <- data.table(X1 = 1, X2 = dtS$X2, X3 = dtS$X3)
X <- as.double(model.matrix(formula(e.T), dtS1))

dprob.deta <- X * exp(-X*coef(e.T)) / (1+ exp(-X*coef(e.T)))^2
pred.logit <- 1 / (1 + exp(-X*coef(e.T)))
iid.logit <-  cbind(dprob.deta) %*% matrix(lava::iid(e.T), nrow = 1)

## GS <- attr(riskRegression:::predictGLM(e.T, newdata = dtS1, average.iid = FALSE), "iid")
## range(iid.logit - GS)

pred.risk <- predict(e.S, newdata = dtS1, times = tau, cause = 1, iid = TRUE)
iid.risk <- pred.risk$absRisk.iid[,1,]

dtS[, pi := pred.logit]
dtS[, r := as.double(pred.risk$absRisk)]

test_that("check point estimate vs. manual calculations", {
    expect_equal(dtS[,mean(r)],e.ateRobust$ate.value["risk.1","Gformula"])

    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.IPCW"])
    expect_equal(dtS[,mean(X1*Y/pi)],e.ateRobust$ate.value["risk.1","IPTW.AIPCW"])

    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(dtS[,mean(X1*Y/pi + r*(1-X1/pi))],e.ateRobust$ate.value["risk.1","AIPTW.AIPCW_knownNuisance"])
})

test_that("check se vs. manual calculations", {

    iid.Gformula1 <- (dtS$r - mean(dtS$r))/NROW(dtS)
    iid.Gformula2 <- colMeans(iid.risk)
    iid.Gformula <- iid.Gformula1 + iid.Gformula2
    expect_equal(sqrt(sum(iid.Gformula^2)),e.ateRobust$ate.se["risk.1","Gformula"])

    iid.IPTW1 <- (dtS[,X1*Y/pi] - dtS[,mean(X1*Y/pi)])/NROW(dtS)
    iid.IPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*dtS$Y/dtS$pi^2) )
    iid.IPTW <- iid.IPTW1 + iid.IPTW2
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.IPCW"])
    expect_equal(sqrt(sum(iid.IPTW^2)),e.ateRobust$ate.se["risk.1","IPTW.AIPCW"])

    iid.AIPTW1 <- (dtS[,X1*Y/pi + r*(1-X1/pi)] - dtS[,mean(X1*Y/pi + r*(1-X1/pi))])/NROW(dtS)
    iid.AIPTW2 <- colMeans( riskRegression::colMultiply_cpp(iid.logit, scale = -dtS$X1*(dtS$Y-dtS$r)/dtS$pi^2) )
    iid.AIPTW3 <- colMeans( riskRegression::colMultiply_cpp(iid.risk, scale = (1-dtS$X1/dtS$pi)) )
    iid.AIPTW <- iid.AIPTW1 + iid.AIPTW2 + iid.AIPTW3
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.IPCW_knownNuisance"])
    expect_equal(sqrt(sum(iid.AIPTW1^2)),e.ateRobust$ate.se["risk.1","AIPTW.AIPCW_knownNuisance"])
})


## ** agreement with ate
set.seed(10)
n <- 1e2
dtS <- sampleData(n,outcome="competing.risks")
## dtS[,min(time),by = event]
tau <- 3 ## tau <- 3:4
e.CSC <- CSC(Hist(time, event) ~ X1 + X2 + X3,
             data = dtS, surv.type = "hazard")

test_that("Agreement ate-ateRobust (competing.risks)",{
    ## NOT POSSIBLE: predictRisk does not recognise the argument product.limit
    ## e.ate <- ate(e.CSC, treatment = "X1", times = 3, data = dtS, cause = 1,
    ## se = TRUE, product.limit = FALSE) 
    e.ateG <- ate(e.CSC, treatment = "X1", times = tau, data = dtS, cause = 1,
                  se = TRUE, verbose = 0)
    e.ateRR <- ate(event = CSC(Hist(time,event) ~ X1 + X2 + X3, data = dtS),
                   treatment = glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")),
                   censor = cph(Surv(time,event==0) ~ X1, data = dtS, x = TRUE, y = TRUE),
                   data = dtS, times = tau, verbose = 0, cause = 1
                   )
    e.ateIP <- ate(event = c("time","event"),
                   treatment = glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")),
                   censor = cph(Surv(time,event==0) ~ X1, data = dtS, x = TRUE, y = TRUE),
                   data = dtS, times = tau, verbose = 0, cause = 1
                   )

    ## data0 <- copy(dtS)
    ## data0[, X1 := factor("0", level = levels(dtS$X1))]
    ## data1 <- copy(dtS)
    ## data1[, X1 := factor("1", level = levels(dtS$X1))]
    ## mean(predict(e.CSC, cause = 1, newdata = data0, times = 3, product.limit = TRUE)$absRisk)
    
    e.ateRobust <- ateRobust(data = dtS, times = tau, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             product.limit = FALSE)
    e.ateRobustPL <- ateRobust(data = dtS, times = tau, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             product.limit = TRUE)

    ## 
    expect_equal(as.double(e.ateRobustPL$ate.value[,"Gformula"]),
                 c(e.ateG$meanRisk[,meanRisk],e.ateG$riskComparison[,diff]),
                 tol = 1e-4)
    expect_equal(as.double(e.ateRobustPL$ate.se[,"Gformula"]),
                 c(e.ateG$meanRisk[,meanRisk.se],e.ateG$riskComparison[,diff.se]),
                 tol = 1e-4)

    expect_equal(as.double(e.ateRobustPL$ate.value[,"IPTW.IPCW"]),
                 c(e.ateIP$meanRisk[,meanRisk],e.ateIP$riskComparison[,diff]),
                 tol = 1e-4)
    ## expect_equal(as.double(e.ateRobustPL$ate.se[,"IPTW.IPCW"]),
                 ## c(e.ateIP$meanRisk[,meanRisk.se],e.ateIP$riskComparison[,diff.se]),
                 ## tol = 1e-4)

    expect_equal(as.double(e.ateRobustPL$ate.value[,"AIPTW.AIPCW_estimatedNuisance"]),
                 c(e.ateRR$meanRisk[,meanRisk],e.ateRR$riskComparison[,diff]),
                 tol = 1e-4)
    ## expect_equal(as.double(e.ateRobustPL$ate.se[,"AIPTW.AIPCW_estimatedNuisance"]),
                 ## c(e.ateRR$meanRisk[,meanRisk.se],e.ateRR$riskComparison[,diff.se]),
                 ## tol = 1e-4)
    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28427, 0.30878, 0.02452), 
                  "IPTW.IPCW" = c(0.25855, 0.5098, 0.25125), 
                  "AIPTW.IPCW_knownNuisance" = c(0.25868, 0.50849, 0.24981), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.25868, 0.50849, 0.24981), 
                  "IPTW.AIPCW" = c(0.25835, 0.50462, 0.24627), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.25848, 0.50331, 0.24483), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.25848, 0.50331, 0.24483)
                  )
    
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5) 
    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.04577, 0.15642, 0.15524), 
                  "IPTW.IPCW" = c(0.04652, 0.16131, 0.16789), 
                  "AIPTW.IPCW_knownNuisance" = c(0.04636, 0.16792, 0.17302), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.04644, 0.15882, 0.16504), 
                  "IPTW.AIPCW" = c(0.04648, 0.1619, 0.16844), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.04632, 0.16799, 0.1731), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0464, 0.15946, 0.16567)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  "IPTW.IPCW" = c(0.25857, 0.5099, 0.25133), 
                  "AIPTW.IPCW_knownNuisance" = c(0.2587, 0.50859, 0.24988), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.2587, 0.50859, 0.24988), 
                  "IPTW.AIPCW" = c(0.25837, 0.50473, 0.24636), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.2585, 0.50341, 0.24491), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.2585, 0.50341, 0.24491)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.04576, 0.15641, 0.15524), 
                  "IPTW.IPCW" = c(0.04653, 0.16134, 0.16792), 
                  "AIPTW.IPCW_knownNuisance" = c(0.04636, 0.16804, 0.17314), 
                  "AIPTW.IPCW_estimatedNuisance" = c(0.04644, 0.15886, 0.16509), 
                  "IPTW.AIPCW" = c(0.04648, 0.16193, 0.16847), 
                  "AIPTW.AIPCW_knownNuisance" = c(0.04632, 0.1681, 0.17322), 
                  "AIPTW.AIPCW_estimatedNuisance" = c(0.0464, 0.15951, 0.16572)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
   
})

## * [ate] Miscellaneous

## ** Pre-computation of iidCox does not affect the results
test_that("[ate] Cox model/G-formula - precompute iid", {
    set.seed(10)
    d <- sampleData(100)

    e.CSC <- CSC(Hist(time,event)~X1+X2, data = d)
    e.CSC <- iidCox(e.CSC)

    GS <- ate(e.CSC,
              treatment = "X1", data = d, times = 1:5, cause = 1,
              verbose = FALSE)
    
    test <- ate(e.CSC,
                treatment = "X1", data = d, times = 1:5, cause = 1,
                verbose = FALSE)

    expect_equal(test$meanRisk, GS$meanRisk)
})

## * [ate] Previous bug
cat("[ate] Previous bug \n")

## ** Results of ate depends on the number of timepoints
test_that("[ate] Cox model/G-formula - one or several timepoints affects the results",{
    d <- sampleData(1000,outcome="survival")
    fit <- coxph(Surv(time,event)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE)
    a1 <- ate(fit,data=d,treatment="X1",time=5,bootci.method="wald", verbose = verbose)
    a2 <- ate(fit,data=d,treatment="X1",time=5:7,bootci.method="wald", verbose = verbose)

    expect_equal(a2$riskComparison[time==5,],
                 a1$riskComparison)
    expect_equal(a2$meanRisk[time==5,],
                 a1$meanRisk[time==5,])
})

## ** Multiple timepoints
## TAG: Monday, Jul 8, 2019 9:09:36 PM (email subject Re: Branche ate for riskRegression ready)
set.seed(10)
n <- 1e2
dtS <- sampleData(n,outcome="competing.risks")
e.CSC <- CSC(Hist(time, event) ~ X1 + X2 + X3,
             data = dtS, surv.type = "hazard")

test_that("ate double robust estimator works with multiple timepoint",{
    ## previous error message
    ## Error in iidTotal[[contrasts[iC]]] + iid.treatment[[1]] : non-numeric argument to binary operator
    e.ateRR <- ate(event = Hist(time,event) ~ X1 + X2 + X3, 
                   treatment = glm(X1 ~ 1, data = dtS, family = binomial(link = "logit")),
                   censor = cph(Surv(time,event==0) ~ X1, data = dtS, x = TRUE, y = TRUE),
                   data = dtS, times = 3:4, verbose = 0, cause = 1
                   )

    e.ateRR <- ate(event = CSC(Hist(time,event) ~ X1 + X2 + X3, data = dtS),
                   treatment = "X1",
                   data = dtS, times = 3:4, verbose = 0, cause = 1
                   )

})
