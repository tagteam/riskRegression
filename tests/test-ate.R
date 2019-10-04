library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(testthat)
library(data.table)
library(ipw)
library(lava)

verbose <- FALSE
calcIterm <- function(factor, iid, indexJump, iid.outsideI){
    n <- length(indexJump)

    if(iid.outsideI){ ## compute iid int factor
        ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
            iIID <- iid[iObs,1,]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else {
                return(iIID*sum(iFactor))
            }
        })
    } else { ## compute int iid * factor
        ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
            iIID <- iid[iObs,1:indexJump[iObs],]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else if(indexJump[iObs]==1){
                return(iIID*iFactor)
            }else{
                return(rowSums(rowMultiply_cpp(t(iIID), iFactor)))
            }
        })
    }

    return(t(do.call(cbind,ls.I)))
}



## * [ate] G-formula only
cat("[ate] G-formula only \n")

## ** data
set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

## ** Compare to fix value / explicit computation 
test_that("[ate] G-formula,survival - compare to fix values / explicit computation",{
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
        ## risk
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

    ## difference in risk
    diffATE <- ATE[["T1"]]-ATE[["T0"]]
    diffATE.iid <- ATE.iid[["T1"]]-ATE.iid[["T0"]]
    diffATE.se <- sqrt(apply(diffATE.iid^2, 2, sum))
    diffATE.lower <- diffATE + qnorm(0.025) * diffATE.se
    diffATE.upper <- diffATE + qnorm(0.975) * diffATE.se

    expect_equal(diffATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.lower])
    expect_equal(diffATE.lower,
                 c(-0.04241346, -0.24364096, -0.26383983), tol = 1e-6)

    expect_equal(diffATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",diff.upper])
    expect_equal(diffATE.upper,
                 c(0.0437218, 0.2512001, 0.2720261), tol = 1e-6)

    ## ratio of the risks
    ratioATE <- ATE[["T1"]]/ATE[["T0"]]
    ratioATE.iid <- rowScale_cpp(ATE.iid[["T1"]],ATE[["T0"]])-rowMultiply_cpp(ATE.iid[["T0"]], ATE[["T1"]]/ATE[["T0"]]^2)
    expect_equal(unname(ATE.iid[["T1"]]),unname(ateFit$iid$T1))

    ratioATE.se <- sqrt(rowSums(ratioATE.iid^2))
    ratioATE.se <- sqrt(colSums(ratioATE.iid^2))
    ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
    ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se
    
    expect_equal(ratioATE.lower,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.lower])
    expect_equal(ratioATE.lower,
                 c(0.2745873, 0.47231036, 0.59979056), tol = 1e-6)
    expect_equal(ratioATE.upper,
                 ateFit$riskComparison[Treatment.A == "T0" & Treatment.B == "T1",ratio.upper])
    expect_equal(ratioATE.upper,
                 c(1.7478076, 1.5439995, 1.41259273))
})

## ** strata argument
test_that("[ate] G-formula,survival - strata argument",{
    ## Cox model
    e.coxph <- coxph(Surv(time, event) ~ X1, data = dtS,
                     x = TRUE, y = TRUE)
    outATE <- ate(e.coxph, data = dtS, treatment = NULL, strata = "X1", time = 1, verbose = verbose)

    outPred <- predictCox(e.coxph,
                          newdata = dtS,
                          times = 1,
                          se = TRUE,
                          keep.newdata = TRUE,
                          type = "survival")
    GS.se <- as.data.table(outPred)[,.SD[1],by = "X1"]
    setkeyv(GS.se, cols = "X1")
    test.se <- outATE$meanRisk
    expect_equal(1-GS.se$survival,test.se$meanRisk)
    expect_equal(GS.se$survival.se,GS.se$survival.se)

    ## Stratified Cox model
    e.coxph <- coxph(Surv(time, event) ~ strata(X1), data = dtS,
                     x = TRUE, y = TRUE)
    outATE <- ate(e.coxph, data = dtS, treatment = NULL, strata = "X1", time = 3, verbose = verbose)

    outPred <- predictCox(e.coxph,
                          newdata = dtS,
                          times = 3,
                          se = TRUE,
                          keep.newdata = TRUE,
                          type = "survival")
    GS.se <- as.data.table(outPred)[,.SD[1],by = "X1"]
    setkeyv(GS.se, cols = "X1")
    test.se <- outATE$meanRisk
    expect_equal(1-GS.se$survival,test.se$meanRisk)
    expect_equal(GS.se$survival.se,GS.se$survival.se)
})    


## * [ate] Logistic regression
cat("[ate] Logistic regression \n")

n <- 100
set.seed(10)
dtB <- sampleData(n, outcome="binary")
## dtB <- data.table(X1 = rbinom(n,size =1,prob=0.5),
                  ## X2 = rbinom(n,size =1,prob=0.5))
## dtB$Y <- rbinom(n, size = 1, prob = expit(dtB$X1+dtB$X2))
## dtB$X1 <- as.factor(dtB$X1)
## dtB$X2 <- as.factor(dtB$X2)

test_that("[ate] logistic regression - compare to lava",{
    fitY <- glm(formula = Y ~ X1 + X2, data=dtB, family = "binomial")
    fitT <- glm(formula = X1 ~ X2, data=dtB, family = "binomial")

    ## G-formula
    e.ate <- ate(fitY, data = dtB, treatment = "X1",
                 times = 5, 
                 se = TRUE, iid = TRUE, B = 0, verbose = verbose)

    e.lava <- estimate(fitY, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        R.X11 <- expit(a + b + c * (data[["X2"]]=="1"))
        R.X10 <- expit(a + c * (data[["X2"]]=="1"))
        list(risk0=R.X10,risk1=R.X11,riskdiff=R.X11-R.X10)},
        average=TRUE)

    expect_equal(unname(e.lava$coef), c(e.ate$meanRisk[,meanRisk], e.ate$riskComparison[,diff]))
    expect_equal(unname(e.lava$vcov[1:2,1:2]), unname(crossprod(do.call(cbind,e.ate$iid))))
    expect_equal(unname(sqrt(diag(e.lava$vcov))), c(e.ate$meanRisk[,meanRisk.se], e.ate$riskComparison[,diff.se]))

    ## AIPTW
    e.ate2 <- ate(fitY,
                  treatment = fitT,
                  data = dtB, 
                  times = 5, 
                  se = TRUE, iid = TRUE, B = 0, verbose = verbose
                  )
    ## ate(fitY,
        ## treatment = fitT,
        ## data = dtB, 
        ## times = 5, 
        ## se = TRUE, iid = TRUE, B = 0, verbose = FALSE, known.nuisance = TRUE
        ## )

    dtB$Y0 <- 0
    dtB$Y1 <- 1

    ## iid ate
    iPredT <- predict(fitT, type = "response")
    dtBC <- rbind(cbind(dtB[,.(Y,X2)], X1 = factor(0, levels = levels(dtB$X1)), X1test = dtB$X1=="0", pi = 1-iPredT),
                  cbind(dtB[,.(Y,X2)], X1 = factor(1, levels = levels(dtB$X1)), X1test = dtB$X1=="1", pi = iPredT))
    dtBC$r <- predict(fitY, newdata = dtBC, type = "response")
    dtBC[, ate := Y*X1test/pi + r*(1-X1test/pi)]
    dtBC[, ate.iid := (ate-mean(ate))/.N, by = "X1"]

    expect_equal(dtBC$ate.iid,
                 do.call(rbind,attr(e.ate2$iid,"ate"))[,1])


    ## iid outcome
    iid.risk <- attr(predictRisk(fitY, newdata = dtBC, iid = TRUE),"iid")
    nuisanceY.iid <- colMultiply_cpp(iid.risk, scale = (1-dtBC$X1test/dtBC$pi))
    dtBC$AnuisanceY.iid <- c(colMeans(nuisanceY.iid[dtBC$X1=="0",]),colMeans(nuisanceY.iid[dtBC$X1=="1",]))

    expect_equal(unname(dtBC$AnuisanceY.iid),
                 unname(do.call(rbind,attr(e.ate2$iid,"outcome"))[,1]))

    ## iid treatment
    iid.pi <- attr(predictRisk(fitT, newdata = dtBC, iid = TRUE),"iid")
    nuisanceT.iid <- colMultiply_cpp(iid.pi, scale = (-1)^(dtBC$X1=="1")*dtBC$X1test*(dtBC$Y-dtBC$r)/dtBC$pi^2)
    dtBC$AnuisanceT.iid <- c(colMeans(nuisanceT.iid[dtBC$X1=="0",]),colMeans(nuisanceT.iid[dtBC$X1=="1",]))
    
    expect_equal(unname(dtBC$AnuisanceT.iid),
                 unname(do.call(rbind,attr(e.ate2$iid,"treatment"))[,1]))
    
    
    ## global
    expect_equal(dtBC[, mean(ate), by = "X1"][[2]],
                 e.ate2$meanRisk[["meanRisk"]])
    ## expect_equal(dtBC[, sqrt(sum(ate.iid^2)), by = "X1"][[2]],
                 ## e.ate2$meanRisk[["meanRisk.se"]])
    expect_equal(dtBC[, sqrt(sum( (ate.iid+AnuisanceY.iid+AnuisanceT.iid)^2 )), by = "X1"][[2]],
                 e.ate2$meanRisk[["meanRisk.se"]])
    ## expect_equal(dtBC[, sqrt(sum( (ate.iid)^2 )), by = "X1"][[2]],
                 ## e.ate2$meanRisk[["meanRisk.se"]])

    
})

## * [ate] Survival case
cat("[ate] Survival case \n")
## ** no censoring - manual computation
## *** Data
n <- 5e1

set.seed(10)
dtS <- sampleData(n,outcome="survival")
tau <- median(dtS$time)

dtS$event <- 1
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)


## *** check
test_that("[ate] no censoring, survival - check vs. manual calculations", {
    e.S <- coxph(Surv(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS,x = TRUE)
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit")) ## dtS$X1

    ## automatic calculation
    e.ate <- list("AIPTW" = ate(data = dtS, times = tau,
                                event = e.S,
                                treatment = e.T,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "AIPTW2" = ate(data = dtS, times = tau,
                                event = e.S,
                                treatment = e.T, method.iid = 2,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW" = ate(data = dtS, times = tau,
                               event = c("time","event"),
                               treatment = e.T,
                               iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW2" = ate(data = dtS, times = tau,
                                event = c("time","event"),
                                treatment = e.T, method.iid = 2,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "G-formula" = ate(data = dtS, times = tau,
                                    event = e.S,
                                    treatment = "X1f",
                                    iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose)
                  )


    expect_equal(e.ate[["IPTW2"]]$meanRisk,e.ate[["IPTW"]]$meanRisk)
    expect_equal(e.ate[["IPTW2"]]$riskComparison,e.ate[["IPTW"]]$riskComparison)
    expect_equal(e.ate[["AIPTW2"]]$meanRisk,e.ate[["AIPTW"]]$meanRisk)
    expect_equal(e.ate[["AIPTW2"]]$riskComparison,e.ate[["AIPTW"]]$riskComparison)
    
    ## manual calculation
    dtCf <- do.call(rbind,lapply(levels(dtS$X1f), function(iT){ ## iT <- "0"
        iDT <- data.table::copy(dtS[,.(event,time,Y,X1,X2,X3,X6)])
        iDT[, X1f := factor(iT,levels(dtS$X1f))]
        iDT[, X1test := iT==as.character(X1)]

        iPred.logit <- predictRisk(e.T, newdata = iDT, iid = TRUE, level = iT)
        iPred.Surv <- predictCox(e.S, newdata = iDT, times = tau, type = "survival", iid = TRUE)

        iDT[, pi := as.double(iPred.logit)]
        iDT[, r := as.double(1-iPred.Surv$survival)]

        iDT[, iid.Gformula := (r - mean(r))/.N]
        iDT[, iid.Gformula := iid.Gformula + colMeans(-iPred.Surv$survival.iid[,1,])]

        iDT[, iid.IPTW := (X1test*Y/pi - mean(X1test*Y/pi))/.N]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y/pi^2))]

        iDT[, iid.AIPTW := (X1test*Y/pi + r*(1-X1test/pi) - mean(X1test*Y/pi + r*(1-X1test/pi)))/.N]
        iDT[, iid.AIPTW := iid.AIPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*(Y-r)/pi^2))]
        iDT[, iid.AIPTW := iid.AIPTW + colMeans(colMultiply_cpp(-iPred.Surv$survival.iid[,1,], scale = 1-X1test/pi))]

        return(iDT)        
    }))
    
    ## check estimate
    expect_equal(dtCf[,mean(r),by="X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.43356267, 0.78109942),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y/pi),by="X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.44100313, 0.74666667),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*(event==1)*(time<=tau)/pi + r*(1-X1test/pi)),by="X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.44729274, 0.76956544),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    ## check influence function
    expect_equal(do.call(rbind,e.ate[["G-formula"]]$iid)[,1],
                 dtCf$iid.Gformula)

    expect_equal(do.call(rbind,e.ate[["IPTW"]]$iid)[,1],
                 dtCf$iid.IPTW)

    expect_equal(do.call(rbind,e.ate[["AIPTW"]]$iid)[,1],
                 dtCf$iid.AIPTW)

    ## check standard error
    expect_equal(sqrt(dtCf[,sum(iid.Gformula^2),by="X1f"][[2]]),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.07227568, 0.09615608),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.IPTW^2),by="X1f"][[2]]),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.07837068, 0.12110765),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.AIPTW^2),by="X1f"][[2]]),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.07701387, 0.10447564),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    ## compare to ipw package
    ww <- ipwpoint(exposure=X1,family="binomial",link="logit",
                   denominator=~X2,
                   data=dtS)$ipw.weights
    expect_equal(ww,
                 dtCf[,1/pi[X1test==1], by = "time"][[2]])
})
    
    

## ** Censoring - manual computation
n <- 5e1

set.seed(10)
dtS <- sampleData(n,outcome="survival")
tau <- median(dtS$time)

dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$C <- (dtS$time<=tau)*(dtS$event!=0)
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))


test_that("[ate] Censoring, survival - check vs. manual calculations", {
    e.S <- coxph(Surv(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, x = TRUE , y = TRUE)
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit"))
    e.C <- coxph(Surv(time, event == 0) ~ X1f + X2, data = dtS, x = TRUE, y = TRUE)
    
    ## automatic calculation
    e.ate <- list("AIPTW" = ate(data = dtS, times = tau,
                                event = e.S,
                                treatment = e.T,
                                censor = e.C,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "AIPTW2" = ate(data = dtS, times = tau,
                                 event = e.S,
                                 treatment = e.T,
                                 censor = e.C, method.iid = 2,
                                 iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW" = ate(data = dtS, times = tau,
                               event = c("time","event"),
                               treatment = e.T,
                               censor = e.C,
                               iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW2" = ate(data = dtS, times = tau,
                                event = c("time","event"),
                                treatment = e.T,
                                censor = e.C, method.iid = 2,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "G-formula" = ate(data = dtS, times = tau,
                                    event = e.S,
                                    treatment = "X1f",
                                    iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose)
                  )


    expect_equal(e.ate[["IPTW2"]]$meanRisk, e.ate[["IPTW"]]$meanRisk)
    expect_equal(e.ate[["IPTW2"]]$riskComparison, e.ate[["IPTW"]]$riskComparison)
    expect_equal(e.ate[["AIPTW2"]]$meanRisk, e.ate[["AIPTW"]]$meanRisk)
    expect_equal(e.ate[["AIPTW2"]]$riskComparison, e.ate[["AIPTW"]]$riskComparison)

    ## ## manual calculation ## ##
    iPred.Cens <- predictCox(e.C, newdata = dtS, times = pmin(tau,dtS$time-1e-10), type = "survival", iid = TRUE)
    iPred.Cens$survivalDiag <- diag(iPred.Cens$survival)
    iPred.Cens$survivalDiag.iid <- t(sapply(1:n, function(iObs){iPred.Cens$survival.iid[iObs,iObs,]}))

    ## augmentation term
    jumpC <- sort(e.C$y[e.C$y[,"status"]==1,"time"])
    indexJump <- prodlim::sindex(jump.time = jumpC, eval.times = pmin(tau,dtS$time))
    pred.Surv_tau <- predictCox(e.S, newdata = dtS, times = tau, type = "survival", iid = TRUE)
    pred.Surv_jump <- predictCox(e.S, newdata = dtS, times = jumpC, type = "survival", iid = TRUE)
    pred.Surv_beforejump <- predictCox(e.S, newdata = dtS, times = jumpC-1e-6, type = "survival", iid = TRUE)
    pred.G_jump <- predictCox(e.C, newdata = dtS, times = jumpC-1e-10, type = "survival", iid = TRUE)
    pred.dN_jump <- do.call(cbind,lapply(jumpC, function(iJ){iJ == dtS$time}))
    pred.dLambda_jump <- predictCox(e.C, newdata = dtS, times = jumpC, type = "hazard", iid = TRUE)
    pred.dM_jump <- pred.dN_jump-pred.dLambda_jump$hazard

    term.SS <- colCenter_cpp(pred.Surv_jump$survival,pred.Surv_tau$survival)/pred.Surv_beforejump$survival
    integrand_jump <- term.SS*pred.dM_jump/pred.G_jump$survival
    augTerm <- sapply(1:n, function(iObs){
        if(indexJump[iObs]==0){
            return(0)
        }else{
            return(sum(integrand_jump[iObs,1:indexJump[iObs]]))
        }
    })
    
    augTermY1.iid <- calcIterm(factor = pred.dM_jump/(pred.Surv_beforejump$survival*pred.G_jump$survival),
                               iid = -pred.Surv_tau$survival.iid,
                               indexJump = indexJump,
                               iid.outsideI = TRUE)

    augTermY2.iid <- calcIterm(factor = -pred.dM_jump/(pred.Surv_beforejump$survival*pred.G_jump$survival),
                               iid = -pred.Surv_jump$survival.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    augTermSurv.iid <- calcIterm(factor = - term.SS * pred.dM_jump/(pred.Surv_beforejump$survival * pred.G_jump$survival),
                                 iid = pred.Surv_beforejump$survival.iid,
                                 indexJump = indexJump,
                                 iid.outsideI = FALSE)
    
    augTermG1.iid <- calcIterm(factor = - term.SS * pred.dM_jump/pred.G_jump$survival^2,
                               iid = pred.G_jump$survival.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    augTermG2.iid <- calcIterm(factor = - term.SS / pred.G_jump$survival,
                               iid = pred.dLambda_jump$hazard.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    dtCf <- do.call(rbind,lapply(levels(dtS$X1f), function(iT){ ## iT <- "0"
        iDT <- data.table::copy(dtS[,.(event,time,Y,C,X1,X2,X3,X6)])
        iDT[, X1f := factor(iT,levels(dtS$X1f))]
        iDT[, X1test := iT==as.character(X1)]

        iPred.logit <- predictRisk(e.T, newdata = iDT, iid = TRUE, level = iT)
        iPred.Surv_tau <- predictCox(e.S, newdata = iDT, times = tau, type = "survival", iid = TRUE)

        iDT[, pi := as.double(iPred.logit)]
        iDT[, r := as.double(1-iPred.Surv_tau$survival)]
        iDT[, G := as.double(iPred.Cens$survivalDiag)]
        iDT[, I := augTerm]

        ## Gformula ##
        iDT[, iid.Gformula := (r - mean(r))/.N]
        iDT[, iid.Gformula := iid.Gformula + colMeans(-iPred.Surv_tau$survival.iid[,1,])]

        ## IPTW
        iDT[, iid.IPTW := (X1test*Y*C/(pi*G) - mean(X1test*Y*C/(pi*G)))/.N]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y*C/(pi^2*G)))]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]

        ## AIPTW,AIPCW ##        
        iDT[, iid.AIPTW.ate := (X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi - mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi))/.N]

        ## outcome
        iDT[, iid.AIPTW.outcome := colMeans(colMultiply_cpp(-iPred.Surv_tau$survival.iid[,1,], (1-X1test/pi)))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + colMeans(colMultiply_cpp(augTermY1.iid,X1test/pi))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + colMeans(colMultiply_cpp(augTermY2.iid,X1test/pi))]

        ## treatment
        iDT[, iid.AIPTW.treatment := colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -(X1test/pi^2)*(Y*C/G - r + I)))]

        ## survival
        iDT[, iid.AIPTW.survival := colMeans(colMultiply_cpp(augTermSurv.iid, X1test/pi))]

        ## censoring
        iDT[, iid.AIPTW.censoring := colMeans(colMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + colMeans(colMultiply_cpp(augTermG1.iid,X1test/pi))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + colMeans(colMultiply_cpp(augTermG2.iid,X1test/pi))]

        ## total
        iDT[, iid.AIPTW := iid.AIPTW.ate + iid.AIPTW.outcome + iid.AIPTW.treatment + iid.AIPTW.survival + iid.AIPTW.censoring]
        return(iDT)        
    }))

    ## check estimate
    expect_equal(dtCf[,mean(r),by="X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.38160373, 0.69088127),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G)),by="X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.37314136, 0.75714469),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi),by="X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.37667566, 0.77294549),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    ## check influence function
    expect_equal(do.call(rbind,e.ate[["G-formula"]]$iid)[,1],
                 dtCf$iid.Gformula)

    expect_equal(do.call(rbind,e.ate[["IPTW"]]$iid)[,1],
                 dtCf$iid.IPTW)

    expect_equal(do.call(rbind,e.ate[["AIPTW"]]$iid)[,1],
                 dtCf$iid.AIPTW)

    ## check standard errors
    expect_equal(dtCf[,sqrt(sum(iid.Gformula^2)), by = "X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.069586, 0.1201558),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.IPTW^2)), by = "X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.07693792, 0.12423590),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.AIPTW^2)), by = "X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.07484163, 0.10880191),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)


})

## * [ate] Competing risk case
cat("[ate] Competing risk case \n")

## ** no censoring - manual computation

## *** Data
n <- 1e2
tau <- 1.5

set.seed(10)
dtS <- sampleData(n, outcome="competing.risks")[event != 0]
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))

## *** check
test_that("[ate] no censoring, competing risks - check vs. manual calculations", {
    e.S <- CSC(Hist(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, surv.type = "survival")
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit")) ## dtS$X1

    ## automatic calculation
    e.ate <- list("AIPTW" = ate(data = dtS, times = tau, cause = 1,
                                event = e.S,
                                treatment = e.T,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "AIPTW2" = ate(data = dtS, times = tau, cause = 1,
                                 event = e.S,
                                 treatment = e.T, method.iid = 2,
                                 iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW" = ate(data = dtS, times = tau, cause = 1,
                               event = c("time","event"),
                               treatment = e.T,
                               iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW2" = ate(data = dtS, times = tau, cause = 1,
                                event = c("time","event"),
                                treatment = e.T, method.iid = 2,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "G-formula" = ate(data = dtS, times = tau, cause = 1,
                                    event = e.S,
                                    treatment = "X1f",
                                    iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose)
                  )

    expect_equal(e.ate[["IPTW2"]]$meanRisk,e.ate[["IPTW"]]$meanRisk)
    expect_equal(e.ate[["IPTW2"]]$riskComparison,e.ate[["IPTW"]]$riskComparison)
    expect_equal(e.ate[["AIPTW2"]]$meanRisk,e.ate[["AIPTW"]]$meanRisk)
    expect_equal(e.ate[["AIPTW2"]]$riskComparison,e.ate[["AIPTW"]]$riskComparison)

    ## manual calculation
    dtCf <- do.call(rbind,lapply(levels(dtS$X1f), function(iT){ ## iT <- "0"
        iDT <- data.table::copy(dtS[,.(event,time,Y,X1,X2,X3,X6)])
        iDT[, X1f := factor(iT,levels(dtS$X1f))]
        iDT[, X1test := iT==as.character(X1)]

        iPred.logit <- predictRisk(e.T, newdata = iDT, iid = TRUE, level = iT)
        iPred.risk <- predict(e.S, newdata = iDT, times = tau, iid = TRUE,
                              cause = 1, product.limit = FALSE)

        iDT[, pi := as.double(iPred.logit)]
        iDT[, r := as.double(iPred.risk$absRisk)]

        iDT[, iid.Gformula := (r - mean(r))/.N]
        iDT[, iid.Gformula := iid.Gformula + colMeans(iPred.risk$absRisk.iid[,1,])]

        iDT[, iid.IPTW := (X1test*Y/pi - mean(X1test*Y/pi))/.N]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y/pi^2))]

        iDT[, iid.AIPTW := (X1test*Y/pi + r*(1-X1test/pi) - mean(X1test*Y/pi + r*(1-X1test/pi)))/.N]
        iDT[, iid.AIPTW := iid.AIPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*(Y-r)/pi^2))]
        iDT[, iid.AIPTW := iid.AIPTW + colMeans(colMultiply_cpp(iPred.risk$absRisk.iid[,1,], scale = 1-X1test/pi))]

        return(iDT)        
    }))

    ## check estimate
    expect_equal(dtCf[,mean(r),by="X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.09315393, 0.37259431),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y/pi),by="X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.08405157, 0.26793249),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*(event==1)*(time<=tau)/pi + r*(1-X1test/pi)),by="X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.07847289, 0.47497984),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    ## check influence function
    expect_equal(do.call(rbind,e.ate[["G-formula"]]$iid)[,1],
                 dtCf$iid.Gformula)

    expect_equal(do.call(rbind,e.ate[["IPTW"]]$iid)[,1],
                 dtCf$iid.IPTW)

    expect_equal(do.call(rbind,e.ate[["AIPTW"]]$iid)[,1],
                 dtCf$iid.AIPTW)

    ## check standard error
    expect_equal(sqrt(dtCf[,sum(iid.Gformula^2),by="X1f"][[2]]),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.0315163, 0.2407176),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.IPTW^2),by="X1f"][[2]]),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.03282852, 0.15183527),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.AIPTW^2),by="X1f"][[2]]),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.03109236, 0.27245093),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)
    
})

## ** Censoring - manual computation
n <- 1e2
tau <- 1.5

set.seed(10)
dtS <- sampleData(n, outcome="competing.risks")
dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$C <- (dtS$time<=tau)*(dtS$event!=0)
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))

test_that("[ate] Censoring, competing risks - check vs. manual calculations", {
    e.S <- CSC(Hist(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, surv.type = "survival")
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit"))
    e.C <- coxph(Surv(time, event == 0) ~ X1f + X2, data = dtS, x = TRUE, y = TRUE)
    
    ## automatic calculation
    e.ate <- list("AIPTW" = ate(data = dtS, times = tau, cause = 1,
                                event = e.S,
                                treatment = e.T,
                                censor = e.C,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "AIPTW2" = ate(data = dtS, times = tau, cause = 1,
                                 event = e.S,
                                 treatment = e.T,
                                 censor = e.C, method.iid = 2,
                                 iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW" = ate(data = dtS, times = tau, cause = 1,
                               event = c("time","event"),
                               treatment = e.T,
                               censor = e.C,
                               iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "IPTW2" = ate(data = dtS, times = tau, cause = 1,
                                event = c("time","event"),
                                treatment = e.T,
                                censor = e.C, method.iid = 2,
                                iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose),
                  "G-formula" = ate(data = dtS, times = tau, cause = 1,
                                    event = e.S,
                                    treatment = "X1f",
                                    iid = TRUE, se = TRUE, product.limit = FALSE, verbose = verbose)
                  )


    expect_equal(e.ate[["IPTW2"]]$meanRisk, e.ate[["IPTW"]]$meanRisk)
    expect_equal(e.ate[["IPTW2"]]$riskComparison, e.ate[["IPTW"]]$riskComparison)
    expect_equal(e.ate[["AIPTW2"]]$meanRisk, e.ate[["AIPTW"]]$meanRisk)
    expect_equal(e.ate[["AIPTW2"]]$riskComparison, e.ate[["AIPTW"]]$riskComparison)

    ## ## manual calculation ## ##
    iPred.Cens <- predictCox(e.C, newdata = dtS, times = pmin(tau,dtS$time-1e-10), type = "survival", iid = TRUE)
    iPred.Cens$survivalDiag <- diag(iPred.Cens$survival)
    iPred.Cens$survivalDiag.iid <- t(sapply(1:n, function(iObs){iPred.Cens$survival.iid[iObs,iObs,]}))

    ## augmentation term
    jumpC <- sort(e.C$y[e.C$y[,"status"]==1,"time"])
    indexJump <- prodlim::sindex(jump.time = jumpC, eval.times = pmin(tau,dtS$time))
    pred.Risk_tau <- predict(e.S, newdata = dtS, times = tau, iid = TRUE, cause = 1)
    pred.Risk_jump <- predict(e.S, newdata = dtS, times = jumpC, iid = TRUE, cause = 1)
    pred.Surv_jump <- predictCox(e.S$models[["OverallSurvival"]], newdata = dtS, times = jumpC-1e-6, type = "survival", iid = TRUE)
    pred.G_jump <- predictCox(e.C, newdata = dtS, times = jumpC-1e-10, type = "survival", iid = TRUE)
    pred.dN_jump <- do.call(cbind,lapply(jumpC, function(iJ){(iJ == dtS$time)*(dtS$event==0)}))
    pred.dLambda_jump <- predictCox(e.C, newdata = dtS, times = jumpC, type = "hazard", iid = TRUE)
    pred.dM_jump <- pred.dN_jump-pred.dLambda_jump$hazard

    term.SS <- -colCenter_cpp(pred.Risk_jump$absRisk,pred.Risk_tau$absRisk)/pred.Surv_jump$survival
    integrand_jump <- term.SS*pred.dM_jump/pred.G_jump$survival
    augTerm <- sapply(1:n, function(iObs){
        if(indexJump[iObs]==0){
            return(0)
        }else{
            return(sum(integrand_jump[iObs,1:indexJump[iObs]]))
        }
    })
    
    augTermY1.iid <- calcIterm(factor = pred.dM_jump/(pred.Surv_jump$survival*pred.G_jump$survival),
                               iid = pred.Risk_tau$absRisk.iid,
                               indexJump = indexJump,
                               iid.outsideI = TRUE)

    augTermY2.iid <- calcIterm(factor = -pred.dM_jump/(pred.G_jump$survival*pred.Surv_jump$survival),
                               iid = pred.Risk_jump$absRisk.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    augTermSurv.iid <- calcIterm(factor = - term.SS * pred.dM_jump/(pred.Surv_jump$survival * pred.G_jump$survival),
                                 iid = pred.Surv_jump$survival.iid,
                                 indexJump = indexJump,
                                 iid.outsideI = FALSE)
    
    augTermG1.iid <- calcIterm(factor = - term.SS * pred.dM_jump/pred.G_jump$survival^2,
                               iid = pred.G_jump$survival.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    augTermG2.iid <- calcIterm(factor = - term.SS / pred.G_jump$survival,
                               iid = pred.dLambda_jump$hazard.iid,
                               indexJump = indexJump,
                               iid.outsideI = FALSE)

    dtCf <- do.call(rbind,lapply(levels(dtS$X1f), function(iT){ ## iT <- "0"
        iDT <- data.table::copy(dtS[,.(event,time,Y,C,X1,X2,X3,X6)])
        iDT[, X1f := factor(iT,levels(dtS$X1f))]
        iDT[, X1test := iT==as.character(X1)]

        iPred.logit <- predictRisk(e.T, newdata = iDT, iid = TRUE, level = iT)
        iPred.risk_tau <- predict(e.S, newdata = iDT, times = tau, iid = TRUE, cause = 1, product.limit = FALSE)

        iDT[, pi := as.double(iPred.logit)]
        iDT[, r := as.double(iPred.risk_tau$absRisk)]
        iDT[, G := as.double(iPred.Cens$survivalDiag)]
        iDT[, I := augTerm]

        ## Gformula ##
        iDT[, iid.Gformula := (r - mean(r))/.N]
        iDT[, iid.Gformula := iid.Gformula + colMeans(iPred.risk_tau$absRisk.iid[,1,])]

        ## IPTW
        iDT[, iid.IPTW := (X1test*Y*C/(pi*G) - mean(X1test*Y*C/(pi*G)))/.N]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y*C/(pi^2*G)))]
        iDT[, iid.IPTW := iid.IPTW + colMeans(colMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]

        ## AIPTW,AIPCW ##        
        iDT[, iid.AIPTW.ate := (X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi - mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi))/.N]

        ## outcome
        iDT[, iid.AIPTW.outcome := colMeans(colMultiply_cpp(iPred.risk_tau$absRisk.iid[,1,], (1-X1test/pi)))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + colMeans(colMultiply_cpp(augTermY1.iid,X1test/pi))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + colMeans(colMultiply_cpp(augTermY2.iid,X1test/pi))]

        ## treatment
        iDT[, iid.AIPTW.treatment := colMeans(colMultiply_cpp(attr(iPred.logit,"iid"), scale = -(X1test/pi^2)*(Y*C/G - r + I)))]

        ## survival
        iDT[, iid.AIPTW.survival := colMeans(colMultiply_cpp(augTermSurv.iid, X1test/pi))]

        ## censoring
        iDT[, iid.AIPTW.censoring := colMeans(colMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + colMeans(colMultiply_cpp(augTermG1.iid,X1test/pi))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + colMeans(colMultiply_cpp(augTermG2.iid,X1test/pi))]

        ## total
        iDT[, iid.AIPTW := iid.AIPTW.ate + iid.AIPTW.outcome + iid.AIPTW.treatment + iid.AIPTW.survival + iid.AIPTW.censoring]
        return(iDT)        
    }))

    ## check estimate
    expect_equal(dtCf[,mean(r),by="X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.07398793, 0.31937921),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G)),by="X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.06717914, 0.19047619),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi),by="X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]])
    expect_equal(c(0.06112147, 0.39625028),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk"]],
                 tol = 1e-6)

    ## check influence function
    expect_equal(do.call(rbind,e.ate[["G-formula"]]$iid)[,1],
                 dtCf$iid.Gformula)

    expect_equal(do.call(rbind,e.ate[["IPTW"]]$iid)[,1],
                 dtCf$iid.IPTW)

    expect_equal(do.call(rbind,e.ate[["AIPTW"]]$iid)[,1],
                 dtCf$iid.AIPTW)

    ## check standard errors
    expect_equal(dtCf[,sqrt(sum(iid.Gformula^2)), by = "X1f"][[2]],
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.02510435, 0.20905295),
                 e.ate[["G-formula"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.IPTW^2)), by = "X1f"][[2]],
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.02644557, 0.12057076),
                 e.ate[["IPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.AIPTW^2)), by = "X1f"][[2]],
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]])
    expect_equal(c(0.02461978, 0.22428431),
                 e.ate[["AIPTW"]]$meanRisk[["meanRisk.se"]],
                 tol = 1e-6)


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
              verbose = verbose)
    
    test <- ate(e.CSC,
                treatment = "X1", data = d, times = 1:5, cause = 1,
                verbose = verbose)

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
                   data = dtS, times = 3:4, verbose = verbose, cause = 1
                   )

    e.ateRR <- ate(event = CSC(Hist(time,event) ~ X1 + X2 + X3, data = dtS),
                   treatment = "X1",
                   data = dtS, times = 3:4, verbose = verbose, cause = 1
                   )

})
