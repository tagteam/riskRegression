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
    dt.ateFit <- as.data.table(ateFit)
    
    ## point estimate
    expect_equal(dt.ateFit[level == "T0" & type == "ate",value],
                 c(0.05842145, 0.46346908, 0.66107744),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0" & type == "ate",lower],
                 c(0.0000000, 0.2756147, 0.4670461),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0" & type == "ate",upper],
                 c(0.1178148, 0.6513235, 0.8551088),
                 tol = 1e-6)

    expect_equal(dt.ateFit[level == "T0.T1" & type == "diffAte",value],
                 c(0.0006541715, 0.0037795576, 0.0040931551),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "diffAte",lower],
                 c(-0.04241346, -0.24364096, -0.26383983),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "diffAte",upper],
                 c(0.0437218, 0.2512001, 0.2720261),
                 tol = 1e-6)

    expect_equal(dt.ateFit[level == "T0.T1" & type == "ratioAte",value],
                 c(1.011197, 1.008155, 1.006192),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "ratioAte",lower],
                 c(0.2745873, 0.4723104, 0.5997906),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "ratioAte",upper],
                 c(1.747808, 1.543999, 1.412593),
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

        expect_equal(pmax(0,ATE.lower), dt.ateFit[level == iT,lower])
        expect_equal(pmin(1,ATE.upper), dt.ateFit[level == iT,upper])
    }

    ## difference in risk
    diffATE <- ATE[["T1"]]-ATE[["T0"]]
    diffATE.iid <- ATE.iid[["T1"]]-ATE.iid[["T0"]]
    diffATE.se <- sqrt(apply(diffATE.iid^2, 2, sum))
    diffATE.lower <- diffATE + qnorm(0.025) * diffATE.se
    diffATE.upper <- diffATE + qnorm(0.975) * diffATE.se

    expect_equal(diffATE.lower, dt.ateFit[level == "T0.T1" & type == "diffAte", lower])
    expect_equal(diffATE.lower, c(-0.04241346, -0.24364096, -0.26383983), tol = 1e-6)

    expect_equal(diffATE.upper, dt.ateFit[level == "T0.T1" & type == "diffAte", upper])
    expect_equal(diffATE.upper, c(0.0437218, 0.2512001, 0.2720261), tol = 1e-6)

    ## ratio of the risks
    ratioATE <- ATE[["T1"]]/ATE[["T0"]]
    ratioATE.iid <- rowScale_cpp(ATE.iid[["T1"]],ATE[["T0"]])-rowMultiply_cpp(ATE.iid[["T0"]], ATE[["T1"]]/ATE[["T0"]]^2)
    expect_equal(unname(ATE.iid[["T1"]]),unname(ateFit$iid$Gformula$T1))

    ratioATE.se <- sqrt(rowSums(ratioATE.iid^2))
    ratioATE.se <- sqrt(colSums(ratioATE.iid^2))
    ratioATE.lower <- ratioATE + qnorm(0.025) * ratioATE.se
    ratioATE.upper <- ratioATE + qnorm(0.975) * ratioATE.se
    
    expect_equal(ratioATE.lower, dt.ateFit[level == "T0.T1" & type == "ratioAte", lower])
    expect_equal(ratioATE.lower, c(0.2745873, 0.47231036, 0.59979056), tol = 1e-6)
    expect_equal(ratioATE.upper, dt.ateFit[level == "T0.T1" & type == "ratioAte", upper])
    expect_equal(ratioATE.upper, c(1.7478076, 1.5439995, 1.41259273))
})

## ** strata argument
test_that("[ate] G-formula,survival - strata argument",{
    ## Cox model
    e.coxph <- coxph(Surv(time, event) ~ X1, data = dtS,
                     x = TRUE, y = TRUE)
    outATE <- ate(e.coxph, data = dtS, treatment = NULL, strata = "X1", time = 1,
                  product.limit = FALSE, verbose = verbose)
    dt.outATE <- as.data.table(outATE)

    outPred <- predictCox(e.coxph,
                          newdata = dtS,
                          times = 1,
                          iid = TRUE,
                          se = TRUE,
                          keep.newdata = TRUE,
                          type = "survival")

    GS.se <- as.data.table(outPred)[,.SD[1],by = "X1"]
    setkeyv(GS.se, cols = "X1")

    expect_equal(1-GS.se$survival, dt.outATE[type=="ate",value])
    expect_equal(GS.se$survival.se, dt.outATE[type=="ate",se])

    ## Stratified Cox model
    e.coxph <- coxph(Surv(time, event) ~ strata(X1), data = dtS,
                     x = TRUE, y = TRUE)
    outATE <- ate(e.coxph, data = dtS, treatment = NULL, strata = "X1", time = 3, verbose = verbose)
    dt.outATE <- as.data.table(outATE)
    
    outPred <- predictCox(e.coxph,
                          newdata = dtS,
                          times = 3,
                          se = TRUE,
                          keep.newdata = TRUE,
                          type = "survival")
    
    GS.se <- as.data.table(outPred)[,.SD[1],by = "X1"]
    setkeyv(GS.se, cols = "X1")

    expect_equal(1-GS.se$survival, dt.outATE[type == "ate",value])
    expect_equal(GS.se$survival.se, dt.outATE[type == "ate",se])
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
    dt.ate <- as.data.table(e.ate)
    
    e.lava <- estimate(fitY, function(p, data){
        a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X21"] ;
        R.X11 <- expit(a + b + c * (data[["X2"]]=="1"))
        R.X10 <- expit(a + c * (data[["X2"]]=="1"))
        list(risk0=R.X10,risk1=R.X11,riskdiff=R.X11-R.X10)},
        average=TRUE)

    expect_equal(unname(e.lava$coef), dt.ate[type %in% c("ate","diffAte"),value])
    expect_equal(unname(e.lava$vcov[1:2,1:2]), unname(crossprod(do.call(cbind,e.ate$iid[["Gformula"]]))))
    expect_equal(unname(sqrt(diag(e.lava$vcov))), c(dt.ate[type %in% c("ate","diffAte"),se]))

    ## AIPTW
    e.ate2 <- ate(fitY,
                  treatment = fitT,
                  data = dtB, 
                  times = 5, 
                  se = TRUE, iid = TRUE, B = 0, verbose = verbose
                  )
    dt.ate2 <- as.data.table(e.ate2)
    
    dtB$Y0 <- 0
    dtB$Y1 <- 1

    ## iid ate
    iPredT <- predict(fitT, type = "response")
    dtBC <- rbind(cbind(dtB[,.(Y,X2)], X1 = factor(0, levels = levels(dtB$X1)), X1test = dtB$X1=="0", pi = 1-iPredT),
                  cbind(dtB[,.(Y,X2)], X1 = factor(1, levels = levels(dtB$X1)), X1test = dtB$X1=="1", pi = iPredT))
    dtBC$r <- predict(fitY, newdata = dtBC, type = "response")
    dtBC[, ate := Y*X1test/pi + r*(1-X1test/pi)]
    dtBC[, ate.iid := (ate-mean(ate))/.N, by = "X1"]

    ## iid outcome
    iid.risk <- attr(predictRisk(fitY, newdata = dtBC, iid = TRUE),"iid")
    nuisanceY.iid <- colMultiply_cpp(iid.risk, scale = (1-dtBC$X1test/dtBC$pi))
    dtBC$AnuisanceY.iid <- c(colMeans(nuisanceY.iid[dtBC$X1=="0",]),colMeans(nuisanceY.iid[dtBC$X1=="1",]))

    ## iid treatment
    iid.pi <- attr(predictRisk(fitT, newdata = dtBC, iid = TRUE),"iid")
    nuisanceT.iid <- colMultiply_cpp(iid.pi, scale = (-1)^(dtBC$X1=="1")*dtBC$X1test*(dtBC$Y-dtBC$r)/dtBC$pi^2)
    dtBC$AnuisanceT.iid <- c(colMeans(nuisanceT.iid[dtBC$X1=="0",]),colMeans(nuisanceT.iid[dtBC$X1=="1",]))
    
    ## global
    expect_equal(dtBC[, mean(ate), by = "X1"][[2]],
                 dt.ate2[type=="ate",value])
    expect_equal(dtBC[, sqrt(sum( (ate.iid+AnuisanceY.iid+AnuisanceT.iid)^2 )), by = "X1"][[2]],
                 dt.ate2[type=="ate",se])
    ## butils:::object2script(dt.ate2)
    ## GS <- data.table("type" = c("ate", "ate", "diffAte", "ratioAte"), 
    ##                  "level" = c("0", "1", "0.1", "0.1"), 
    ##                  "time" = c(5, 5, 5, 5), 
    ##                  "value" = c(0.4167808, 0.7341667, 0.3173858, 1.7615174), 
    ##                  "se" = c(0.05227194, 0.13439533, 0.14429419, 0.39122437), 
    ##                  "lower" = c(0.31432969, 0.47075666, 0.03457442, 0.99473171), 
    ##                  "upper" = c(0.5192320, 0.9975767, 0.6001973, 2.5283031), 
    ##                  "estimator" = c("AIPTW", "AIPTW", "AIPTW", "AIPTW"))
    
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
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, se = TRUE, product.limit = FALSE
                  )
    dt.ate1 <- as.data.table(e.ate1)
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, se = TRUE, band = TRUE, product.limit = FALSE
                  )
    dt.ate2 <- as.data.table(e.ate2)

    expect_equal(dt.ate1, dt.ate2[,.SD,.SDcols = names(dt.ate1)], tol = 1e-8)
    
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
                 dt.ate1[estimator == "Gformula" & type == "ate",value])
    expect_equal(c(0.43356267, 0.78109942),
                 dt.ate1[estimator == "Gformula" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y/pi),by="X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",value])
    expect_equal(c(0.44100313, 0.74666667),
                 dt.ate1[estimator == "IPTW" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*(event==1)*(time<=tau)/pi + r*(1-X1test/pi)),by="X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",value])
    expect_equal(c(0.44729274, 0.76956544),
                 dt.ate1[estimator == "AIPTW" & type == "ate",value],
                 tol = 1e-6)

    ## check influence function
    expect_equal(unname(do.call(rbind,e.ate1$iid[["Gformula"]])[,1]),
                 unname(dtCf$iid.Gformula))

    expect_equal(unname(do.call(rbind,e.ate1$iid[["IPTW"]])[,1]),
                 unname(dtCf$iid.IPTW))

    expect_equal(unname(do.call(rbind,e.ate1$iid[["AIPTW"]])[,1]),
                 unname(dtCf$iid.AIPTW))

    ## check standard error
    expect_equal(sqrt(dtCf[,sum(iid.Gformula^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "Gformula" & type == "ate",se])
    expect_equal(c(0.07227568, 0.09615608),
                 dt.ate1[estimator == "Gformula" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.IPTW^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "IPTW" & type == "ate",se])
    expect_equal(c(0.07837068, 0.12110765),
                 dt.ate1[estimator == "IPTW" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.AIPTW^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se])
    expect_equal(c(0.07701387, 0.10447564),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se],
                 tol = 1e-6)

    ##  check transformation
    GS <- data.table("value" = c(0.3475368, 0.3056635, 0.3222727), 
                     "se" = c(0.1028439, 0.1420700, 0.1109526), 
                     "lower" = c(0.14596632, 0.02721154, 0.10480955), 
                     "upper" = c(0.5491072, 0.5841155, 0.5397359), 
                     "p.value" = c(0.000726804, 0.031436739, 0.003677258), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(1.801584, 1.693110, 1.720496), 
                     "se" = c(0.3194896, 0.4006286, 0.3214094), 
                     "lower" = c(1.175396, 0.907892, 1.090545), 
                     "upper" = c(2.427772, 2.478327, 2.350447), 
                     "p.value" = c(0.01210904, 0.08362040, 0.02498223), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,estimator)],
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
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  censor = e.C,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE
                  )
    dt.ate1 <- as.data.table(e.ate1)
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  censor = e.C,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE
                  )
    dt.ate2 <- as.data.table(e.ate2)

    expect_equal(dt.ate1, dt.ate2, tol = 1e-8)

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
                 dt.ate1[estimator == "Gformula" & type == "ate",value])
    expect_equal(c(0.38160373, 0.69088127),
                 dt.ate1[estimator == "Gformula" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G)),by="X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",value])
    expect_equal(c(0.37314136, 0.75714469),
                 dt.ate1[estimator == "IPTW" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi),by="X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",value])
    expect_equal(c(0.37667566, 0.77294549),
                 dt.ate1[estimator == "AIPTW" & type == "ate",value],
                 tol = 1e-6)

    ## check influence function
    expect_equal(unname(do.call(rbind,e.ate1$iid[["Gformula"]])[,1]),
                 dtCf$iid.Gformula)

    expect_equal(unname(do.call(rbind,e.ate1$iid[["IPTW"]])[,1]),
                 dtCf$iid.IPTW)

    expect_equal(unname(do.call(rbind,e.ate1$iid[["AIPTW"]])[,1]),
                 dtCf$iid.AIPTW)

    ## check standard errors
    expect_equal(dtCf[,sqrt(sum(iid.Gformula^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "Gformula" & type == "ate",se])
    expect_equal(c(0.069586, 0.1201558),
                 dt.ate1[estimator == "Gformula" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.IPTW^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",se])
    expect_equal(c(0.07693792, 0.12423590),
                 dt.ate1[estimator == "IPTW" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.AIPTW^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",se])
    expect_equal(c(0.07484163, 0.10880191),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se],
                 tol = 1e-6)

    ##  check transformation
    GS <- data.table("value" = c(0.3092775, 0.3840033, 0.3962698), 
                     "se" = c(0.1227820, 0.1425698, 0.1121351), 
                     "lower" = c(0.0686292, 0.1045717, 0.1764890), 
                     "upper" = c(0.5499259, 0.6634350, 0.6160507), 
                     "p.value" = c(0.0117716889, 0.0070718645, 0.0004095391), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(1.810468, 2.029110, 2.052019), 
                     "se" = c(0.3948003, 0.5204967, 0.4234462), 
                     "lower" = c(1.036673, 1.008955, 1.222079), 
                     "upper" = c(2.584262, 3.049264, 2.881958), 
                     "p.value" = c(0.04008663, 0.04802263, 0.01297622), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS,dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,estimator)],
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
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose, cause = 1,
                  iid = TRUE, product.limit = FALSE
                  )
    dt.ate1 <- as.data.table(e.ate1)
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose, cause = 1,
                  iid = TRUE, product.limit = FALSE
                  )
    dt.ate2 <- as.data.table(e.ate2)

    expect_equal(dt.ate1,  dt.ate2)

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
                 dt.ate1[estimator == "Gformula" & type == "ate",value])
    expect_equal(c(0.09315393, 0.37259431),
                 dt.ate1[estimator == "Gformula" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y/pi),by="X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",value])
    expect_equal(c(0.08405157, 0.26793249),
                 dt.ate1[estimator == "IPTW" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*(event==1)*(time<=tau)/pi + r*(1-X1test/pi)),by="X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",value])
    expect_equal(c(0.07847289, 0.47497984),
                 dt.ate1[estimator == "AIPTW" & type == "ate",value],
                 tol = 1e-6)

    ## check influence function
    expect_equal(unname(do.call(rbind,e.ate1$iid[["Gformula"]])[,1]),
                 unname(dtCf$iid.Gformula))

    expect_equal(unname(do.call(rbind,e.ate1$iid[["IPTW"]])[,1]),
                 unname(dtCf$iid.IPTW))

    expect_equal(unname(do.call(rbind,e.ate1$iid[["AIPTW"]])[,1]),
                 unname(dtCf$iid.AIPTW))

    ## check standard error
    expect_equal(sqrt(dtCf[,sum(iid.Gformula^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "Gformula" & type == "ate",se])
    expect_equal(c(0.0315163, 0.2407176),
                 dt.ate1[estimator == "Gformula" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.IPTW^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "IPTW" & type == "ate",se])
    expect_equal(c(0.03282852, 0.15183527),
                 dt.ate1[estimator == "IPTW" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(sqrt(dtCf[,sum(iid.AIPTW^2),by="X1f"][[2]]),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se])
    expect_equal(c(0.03109236, 0.27245093),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se],
                 tol = 1e-6)
    
    ##  check transformation
    GS <- data.table("value" = c(0.2794404, 0.1838809, 0.3965069), 
                     "se" = c(0.2272017, 0.1549114, 0.2688568), 
                     "lower" = c(-0.1658668, -0.1197399, -0.1304426), 
                     "upper" = c(0.7247475, 0.4875017, 0.9234565), 
                     "p.value" = c(0.2187263, 0.2352249, 0.1402693), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(3.999770, 3.187715, 6.052789), 
                     "se" = c(2.266232, 2.180114, 3.865633), 
                     "lower" = c(0, 0, 0), 
                     "upper" = c( 8.441504,  7.460660, 13.629289), 
                     "p.value" = c(0.1856088, 0.3156261, 0.1911770), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

})

## ** Censoring - manual computation, surv.type="survival"
n <- 1e2

set.seed(10)
dtS <- sampleData(n, outcome="competing.risks")
tau <- median(dtS$time)

dtS$Y <- (dtS$time<=tau)*(dtS$event==1)
dtS$C <- (dtS$time<=tau)*(dtS$event!=0)
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))

test_that("[ate] Censoring, competing risks (surv.type=\"survival\") - check vs. manual calculations", {
    e.S <- CSC(Hist(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, surv.type = "survival")
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit"))
    e.C <- coxph(Surv(time, event == 0) ~ X1f + X2, data = dtS, x = TRUE, y = TRUE)
    
    ## automatic calculation
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  censor = e.C,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE, cause = 1
                  )
    dt.ate1 <- as.data.table(e.ate1)
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  censor = e.C,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE, cause = 1
                  )
    dt.ate2 <- as.data.table(e.ate2)

    expect_equal(dt.ate1, dt.ate2, tol = 1e-8)
    
    ## ## manual calculation ## ##
    iPred.Cens <- predictCox(e.C, newdata = dtS, times = pmin(tau,dtS$time-1e-10), type = "survival", iid = TRUE)
    iPred.Cens$survivalDiag <- diag(iPred.Cens$survival)
    iPred.Cens$survivalDiag.iid <- t(sapply(1:n, function(iObs){iPred.Cens$survival.iid[iObs,iObs,]}))

    ## augmentation term
    jumpC <- sort(e.C$y[e.C$y[,"status"]==1,"time"])
    indexJump <- prodlim::sindex(jump.time = jumpC, eval.times = pmin(tau,dtS$time))
    pred.Risk_tau <- predict(e.S, newdata = dtS, times = tau, iid = TRUE, cause = 1, product.limit = FALSE)
    pred.Risk_jump <- predict(e.S, newdata = dtS, times = jumpC, iid = TRUE, cause = 1, product.limit = FALSE)
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
                 dt.ate1[estimator == "Gformula" & type == "ate",value])
    expect_equal(c(0.32111469, 0.99623763),
                 dt.ate1[estimator == "Gformula" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G)),by="X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",value])
    expect_equal(c(0.31871273, 0.49590715),
                 dt.ate1[estimator == "IPTW" & type == "ate",value],
                 tol = 1e-6)

    expect_equal(dtCf[,mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi),by="X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",value])
    expect_equal(c(0.29549588, 0.9961723),
                 dt.ate1[estimator == "AIPTW" & type == "ate",value],
                 tol = 1e-6)
    
    ## check influence function
    expect_equal(unname(do.call(rbind,e.ate1$iid[["Gformula"]])[,1]),
                 dtCf$iid.Gformula)

    expect_equal(unname(do.call(rbind,e.ate1$iid[["IPTW"]])[,1]),
                 dtCf$iid.IPTW)
    
    expect_equal(unname(do.call(rbind,e.ate1$iid[["AIPTW"]])[,1]),
                 dtCf$iid.AIPTW)

    ## check standard errors
    expect_equal(dtCf[,sqrt(sum(iid.Gformula^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "Gformula" & type == "ate",se])
    expect_equal(c(0.05007538, 0.42597332),
                 dt.ate1[estimator == "Gformula" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.IPTW^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "IPTW" & type == "ate",se])
    expect_equal(c(0.04961228, 0.15883446),
                 dt.ate1[estimator == "IPTW" & type == "ate",se],
                 tol = 1e-6)

    expect_equal(dtCf[,sqrt(sum(iid.AIPTW^2)), by = "X1f"][[2]],
                 dt.ate1[estimator == "AIPTW" & type == "ate",se])
    expect_equal(c(0.04803356, 0.316899),
                 dt.ate1[estimator == "AIPTW" & type == "ate",se],
                 tol = 1e-6)
    
    ##  check transformation
    GS <- data.table("value" = c(0.6751229, 0.1771944, 0.7006764), 
                     "se" = c(0.4352720, 0.1655279, 0.3273409), 
                     "lower" = c(-0.17799448, -0.14723439,  0.05910013), 
                     "upper" = c(1.0000000, 0.5016232, 1.0000000), 
                     "p.value" = c(0.12089283, 0.28440314, 0.03231356), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(3.102436, 1.555969, 3.371188), 
                     "se" = c(1.4694536, 0.5500765, 1.2732042), 
                     "lower" = c(0.2223595, 0.4778388, 0.8757540), 
                     "upper" = c(5.982512, 2.634099, 5.866623), 
                     "p.value" = c(0.15249898, 0.31215422, 0.06254973), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,estimator)],
                 tol = 1e-6)

})


## ** Censoring - manual computation, surv.type="hazard"
n <- 1e2

set.seed(10)
dtS <- sampleData(n, outcome="competing.risks")
tau <- median(dtS$time)

dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))

test_that("[ate] Censoring, competing risks (surv.type=\"hazard\") - check vs. manual calculations", {
    e.S <- CSC(Hist(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, surv.type = "hazard")
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit"))
    e.C <- coxph(Surv(time, event == 0) ~ X1f + X2, data = dtS, x = TRUE, y = TRUE)
    
    ## automatic calculation
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  censor = e.C,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE, cause = 1
                  )
    dt.ate1 <- as.data.table(e.ate1)
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  censor = e.C,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  iid = TRUE, product.limit = FALSE, cause = 1
                  )
    dt.ate2 <- as.data.table(e.ate2)

    expect_equal(dt.ate1, dt.ate2, tol = 1e-8)
})

## * [ate] Landmark analysis (time varying covariates)
cat("[ate] Landmark analysis \n")
data(veteran, package = "survival")

fit <- coxph(Surv(time, status) ~ celltype+karno + age + trt, veteran)
vet2 <- survSplit(Surv(time, status) ~., veteran,
                  cut=c(60, 120), episode ="timegroup")
fitTD <- coxph(Surv(tstart, time, status) ~ celltype+karno + age + trt,
               data= vet2,x=1)

test_that("[ate] landmark analyses", {
    resVet <- ate(fitTD,formula=Hist(entry=tstart,time=time,event=status)~1,
                  data = vet2, treatment = "celltype", contrasts = NULL,
                  times=5,verbose=verbose,
                  landmark = c(0,30,60,90), cause = 1, se = FALSE)
    dt.resVet <- as.data.table(resVet)
    
    GS <- c(0.01773158, 0.02092285, 0.01489588, 0.02981096, 0.04098041, 0.04821725, 0.03463358, 0.06846231, 0.05589160, 0.06564467, 0.04741838, 0.09298832, 0.02633366, 0.03103968, 0.02217127, 0.04416996)
    expect_equal(dt.resVet[type=="ate",value], GS, tol = 1e-6)

    GS <- c(0.02324883, 0.0272944, 0.01973769, 0.03865134, 0.03816002, 0.04472182, 0.0325225, 0.06317735, 0.00860208, 0.01011683, 0.00727539, 0.014359, 0.01491119, 0.01742742, 0.0127848, 0.02452601, -0.01464675, -0.01717757, -0.01246231, -0.02429234, -0.02955794, -0.03460499, -0.02524711, -0.04881836)
    expect_equal(dt.resVet[type=="diffAte",value], GS, tol = 1e-6)

    GS <- c(2.31115372, 2.30452621, 2.32504339, 2.29654787, 3.15209352, 3.13746333, 3.18332093, 3.11926569, 1.48512762, 1.4835302, 1.488416, 1.48166834, 1.36386147, 1.3614353, 1.36914474, 1.3582411, 0.64259145, 0.64374629, 0.64016698, 0.64517198, 0.47115595, 0.47284384, 0.46756706, 0.47500549)
    expect_equal(dt.resVet[type=="ratioAte",value], GS, tol = 1e-6)

})
## * [ate] Miscellaneous

## ** Pre-computation of iidCox does not affect the results
test_that("[ate] Cox model/G-formula - precompute iid", {
    set.seed(10)
    d <- sampleData(100)

    e.CSC <- CSC(Hist(time,event)~X1+X2, data = d)
    
    GS <- ate(e.CSC,
              treatment = "X1", data = d, times = 1:5, cause = 1,
              verbose = verbose)

    e.CSC <- iidCox(e.CSC)

    test <- ate(e.CSC,
                treatment = "X1", data = d, times = 1:5, cause = 1,
                verbose = verbose)

    expect_equal(as.data.table(test), as.data.table(GS))
})
## ** Bootstrap
if(FALSE){
n <- 250

set.seed(10)
dtS <- sampleData(n,outcome="survival")
tau <- median(dtS$time)
dtS$X1f <- dtS$X1
dtS$X1 <- as.numeric(as.character(dtS$X1))

test_that("[ate] check that bootstrap returns a result", {
    e.S <- coxph(Surv(time, event) ~ X1f + strata(X2) + X3*X6, data = dtS, x = TRUE , y = TRUE)
    e.T <- glm(X1f ~ X2, data = dtS, family = binomial(link = "logit"))
    e.C <- coxph(Surv(time, event == 0) ~ X1f + X2, data = dtS, x = TRUE, y = TRUE)

    e.GS <- ate(data = dtS, times = tau,
                event = e.S,
                treatment = e.T,
                censor = e.C,
                estimator = c("Gformula","IPTW","AIPTW"),
                verbose = verbose,
                se = TRUE, iid = FALSE, product.limit = FALSE,
                )
    
    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  censor = e.C,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  se = TRUE, iid = FALSE, product.limit = FALSE,
                  B = 1000, cpus = 1
                  )

    rbind(e.GS$meanRisk,e.ate1$meanRisk[,.SD,.SDcols = names(e.GS$meanRisk)])
    rbind(e.GS$riskComparison,e.ate1$riskComparison[,.SD,.SDcols = names(e.GS$riskComparison)])
    dt.ate1 <- as.data.table(e.ate1)
})
}
## * [ate] Previous bug
cat("[ate] Previous bug \n")

## ** Results of ate depends on the number of timepoints
test_that("[ate] Cox model/G-formula - one or several timepoints affects the results",{
    d <- sampleData(1000,outcome="survival")
    fit <- coxph(Surv(time,event)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE)
    a1 <- ate(fit,data=d,treatment="X1",time=5, verbose = verbose)
    a2 <- ate(fit,data=d,treatment="X1",time=5:7, verbose = verbose)

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
