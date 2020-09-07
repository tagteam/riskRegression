## minor change compare to the hard coded version can be due to a previous bug in iidCox (the value have not been updated)
## (the calculation of the iid was stopped on jump to early when using tau.max)
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
            iIID <- iid[,1,iObs]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else {
                return(iIID*sum(iFactor))
            }
        })
    } else { ## compute int iid * factor
        ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
            iIID <- iid[,1:indexJump[iObs],iObs]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else if(indexJump[iObs]==1){
                return(iIID*iFactor)
            }else{
                return(rowSums(rowMultiply_cpp(iIID, iFactor)))
            }
        })
    }
    return(do.call(cbind,ls.I))
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
                  times = seqTime, band = FALSE, verbose = verbose)
    dt.ateFit <- as.data.table(confint(ateFit, band = TRUE, p.value = FALSE))

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

    expect_equal(dt.ateFit[level == "T0.T1" & type == "diffAte", lowerBand],
                 c(-0.04230655, -0.24302677, -0.26317472),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "diffAte", upperBand],
                 c(0.04361489, 0.25058589, 0.27136103),
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

    expect_equal(dt.ateFit[level == "T0.T1" & type == "ratioAte", lowerBand],
                 c(0.26976028, 0.46879896, 0.5971274),
                 tol = 1e-6)
    expect_equal(dt.ateFit[level == "T0.T1" & type == "ratioAte", upperBand],
                 c(1.75263463, 1.5475109, 1.41525588),
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
        
        ATE.iid_term1 <- apply(-resPred$survival.iid,1,rowMeans)
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
                  product.limit = FALSE, verbose = verbose, band = TRUE)
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
    outATE <- ate(e.coxph, data = dtS, treatment = NULL, strata = "X1", time = 3,
                  verbose = verbose, band = TRUE)
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
    e0.ate <- ate(fitY, data = dtB, treatment = "X1",
                 se = TRUE, iid = TRUE, B = 0, verbose = verbose)
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
    nuisanceY.iid <- rowMultiply_cpp(iid.risk, scale = (1-dtBC$X1test/dtBC$pi))
    dtBC$AnuisanceY.iid <- c(rowMeans(nuisanceY.iid[,dtBC$X1=="0"]),rowMeans(nuisanceY.iid[,dtBC$X1=="1"]))

    ## iid treatment
    iid.pi <- attr(predictRisk(fitT, newdata = dtBC, iid = TRUE),"iid")
    nuisanceT.iid <- rowMultiply_cpp(iid.pi, scale = (-1)^(dtBC$X1=="1")*dtBC$X1test*(dtBC$Y-dtBC$r)/dtBC$pi^2)
    dtBC$AnuisanceT.iid <- c(rowMeans(nuisanceT.iid[,dtBC$X1=="0"]),rowMeans(nuisanceT.iid[,dtBC$X1=="1"]))
    
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
                  band = FALSE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate1 <- as.data.table(confint(e.ate1, band = TRUE))
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  band = FALSE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate2 <- as.data.table(confint(e.ate2, band = TRUE))

    expect_equal(dt.ate1[,.SD,.SDcols = names(dt.ate2)], dt.ate2, tol = 1e-8)
    
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
        iDT[, iid.Gformula := iid.Gformula + rowMeans(-iPred.Surv$survival.iid[,1,])]

        iDT[, iid.IPTW := (X1test*Y/pi - mean(X1test*Y/pi))/.N]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y/pi^2))]

        iDT[, iid.AIPTW := (X1test*Y/pi + r*(1-X1test/pi) - mean(X1test*Y/pi + r*(1-X1test/pi)))/.N]
        iDT[, iid.AIPTW := iid.AIPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*(Y-r)/pi^2))]
        iDT[, iid.AIPTW := iid.AIPTW + rowMeans(rowMultiply_cpp(-iPred.Surv$survival.iid[,1,], scale = 1-X1test/pi))]

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
                     "p.value" = c(0.00072680, 0.03143674, 0.00367726), 
                     "quantileBand" = c(1.962359, 1.973331, 1.915003), 
                     "adj.p.value" = c(0.0007, 0.0296, 0.0056), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(1.801584, 1.693110, 1.720496), 
                     "se" = c(0.3194896, 0.4006286, 0.3214094), 
                     "lower" = c(1.175396, 0.907892, 1.090545), 
                     "upper" = c(2.427772, 2.478327, 2.350447), 
                     "p.value" = c(0.01210904, 0.08362040, 0.02498223), 
                     "quantileBand" = c(1.971454, 1.969537, 1.952972), 
                     "adj.p.value" = c(0.0118, 0.0832, 0.0246), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
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
                  band = FALSE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate1 <- as.data.table(confint(e.ate1, band = TRUE, p.value = TRUE))
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  censor = e.C,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  band = FALSE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate2 <- as.data.table(confint(e.ate2, band = TRUE, p.value = TRUE))

    expect_equal(dt.ate1[,.SD,.SDcols = names(dt.ate2)], dt.ate2, tol = 1e-8)

    ## ## manual calculation ## ##
    iPred.Cens <- predictCox(e.C, newdata = dtS, times = pmin(tau,dtS$time-1e-10), type = "survival", iid = TRUE)
    iPred.Cens$survivalDiag <- diag(iPred.Cens$survival)
    iPred.Cens$survivalDiag.iid <- sapply(1:n, function(iObs){iPred.Cens$survival.iid[,iObs,iObs]})

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
        iDT[, iid.Gformula := iid.Gformula + rowMeans(-iPred.Surv_tau$survival.iid[,1,])]

        ## IPTW
        iDT[, iid.IPTW := (X1test*Y*C/(pi*G) - mean(X1test*Y*C/(pi*G)))/.N]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y*C/(pi^2*G)))]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]

        ## AIPTW,AIPCW ##        
        iDT[, iid.AIPTW.ate := (X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi - mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi))/.N]

        ## outcome
        iDT[, iid.AIPTW.outcome := rowMeans(rowMultiply_cpp(-iPred.Surv_tau$survival.iid[,1,], (1-X1test/pi)))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + rowMeans(rowMultiply_cpp(augTermY1.iid,X1test/pi))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + rowMeans(rowMultiply_cpp(augTermY2.iid,X1test/pi))]

        ## treatment
        iDT[, iid.AIPTW.treatment := rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -(X1test/pi^2)*(Y*C/G - r + I)))]

        ## survival
        iDT[, iid.AIPTW.survival := rowMeans(rowMultiply_cpp(augTermSurv.iid, X1test/pi))]

        ## censoring
        iDT[, iid.AIPTW.censoring := rowMeans(rowMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + rowMeans(rowMultiply_cpp(augTermG1.iid,X1test/pi))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + rowMeans(rowMultiply_cpp(augTermG2.iid,X1test/pi))]

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
                     "se" = c(0.1227820, 0.1425698, 0.1121352), 
                     "lower" = c(0.0686292, 0.1045717, 0.1764890), 
                     "upper" = c(0.5499259, 0.6634350, 0.6160507), 
                     "p.value" = c(0.01177169, 0.00707186, 0.00040954), 
                     "quantileBand" = c(1.965073, 1.948196, 1.931335), 
                     "adj.p.value" = c(0.0122, 0.0079, 0.0007), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(1.810468, 2.029110, 2.052019), 
                     "se" = c(0.3948003, 0.5204967, 0.4234462), 
                     "lower" = c(1.036673, 1.008955, 1.222079), 
                     "upper" = c(2.584262, 3.049264, 2.881958), 
                     "p.value" = c(0.04008663, 0.04802263, 0.01297622), 
                     "quantileBand" = c(1.961192, 1.959446, 1.952883), 
                     "adj.p.value" = c(0.0392, 0.0459, 0.0129), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS,dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
                 tol = 1e-6)

})

## * [ate] IPCW logistic regression
set.seed(10)
n <- 500
tau <- 1:5
d <- sampleData(n, outcome = "competing.risks")
d$Y <- (d$event == 1)*(d$time <= tau[3])
d0 <- d[event!=0] ## remove censoring

test_that("[ate] IPCW LR - no censoring", {
    e0.wglm <- wglm(regressor.event = ~ X1+X2+X6, formula.censor = Surv(time,event==0) ~ X1,
                    times = tau, data = d0)
    e0.glm <- glm(Y ~ X1+X2+X6, family = binomial, data = d0)

    ## NOTE difference between wglm and glm is that lava:::information.glm uses numerical derivatives
    ## while information.wglm uses explicit formula
    ## this lead to small difference in the influence function and therefore in the standard error
    eATE.wglm <- ate(e0.wglm, data = d0, treatment = "X1", times = tau, band = FALSE, verbose = FALSE)
    eATE.glm <- ate(e0.glm, data = d0, treatment = "X1", times = tau[3], verbose = FALSE)
    expect_equal(eATE.glm$meanRisk, eATE.wglm$meanRisk[time==tau[3]], tol = 1e-4) 

    eATE2.wglm <- ate(e0.wglm, data = d0, treatment = X1~X2, times = tau, band = FALSE, verbose = FALSE)
    eATE2.glm <- ate(e0.glm, data = d0, treatment = X1~X2, times = tau[3], verbose = FALSE)
    expect_equal(eATE2.glm$meanRisk, eATE2.wglm$meanRisk[time==tau[3]], tol = 1e-4)
})

## ** AIPTW
## *** data
d <- data.frame("Ytrue" = c(0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
           "status" = c(2, 2, 0, 1, 1, 0, 1, 0, 0, 2, 2, 0, 0, 1, 0, 1, 2, 0, 2, 1, 2, 1, 2, 0, 1, 1, 0, 1, 2, 0, 1, 2, 1, 0, 0, 1, 2, 1, 1, 0, 1, 2, 2, 1, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 1, 0, 1, 1, 2, 0, 0, 1, 1, 1, 1, 1, 1, 0, 2, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 0, 1, 1, 1, 1, 2, 1, 1, 1, 0, 2, 2, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 2, 0, 0, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 0, 1, 1, 2, 0, 2, 2, 1, 0, 1, 2, 0, 1, 2, 0, 2, 2, 0, 2, 0, 1, 0, 1, 0, 1, 2, 2, 0, 2, 2, 2, 0, 2, 2, 1, 1, 0, 2, 2, 1, 1, 1, 1, 2, 0, 1, 0, 2, 1, 0, 0, 2, 2, 2, 1, 2, 0, 1, 2, 1, 1, 2, 0, 0, 2, 1, 2, 0, 0, 2, 2, 2, 1, 1, 0, 0, 0, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1), 
           "time" = c(0.1029, 0.133, 0.1403, 0.1659, 0.171, 0.1868, 0.1973, 0.2055, 0.2083, 0.2182, 0.2598, 0.3054, 0.3277, 0.3295, 0.3558, 0.3581, 0.3697, 0.3952, 0.4796, 0.4884, 0.5139, 0.6098, 0.619, 0.6255, 0.6374, 0.6517, 0.6687, 0.6856, 0.7404, 0.7557, 0.7625, 0.7757, 0.831, 0.8844, 0.9163, 0.9347, 0.9645, 1.0534, 1.078, 1.0848, 1.1464, 1.1603, 1.1778, 1.179, 1.1942, 1.1948, 1.2169, 1.2453, 1.2701, 1.2899, 1.3211, 1.3882, 1.413, 1.4392, 1.4407, 1.4407, 1.441, 1.451, 1.5042, 1.5496, 1.5595, 1.5974, 1.8156, 1.8209, 1.8228, 1.8495, 1.8504, 1.8515, 1.8887, 1.9146, 1.915, 2.0289, 2.0512, 2.0874, 2.1309, 2.1318, 2.1405, 2.1658, 2.1674, 2.1797, 2.1944, 2.2099, 2.2665, 2.2956, 2.4062, 2.444, 2.4472, 2.4835, 2.5011, 2.5401, 2.5416, 2.5741, 2.5825, 2.6255, 2.6335, 2.6452, 2.649, 2.6695, 2.7385, 2.7989, 2.8758, 2.8894, 2.9074, 2.9453, 2.996, 3.0597, 3.0712, 3.0835, 3.0847, 3.1286, 3.1626, 3.163, 3.1785, 3.182, 3.2286, 3.3285, 3.3329, 3.4164, 3.4316, 3.5536, 3.5564, 3.5627, 3.6225, 3.6551, 3.668, 3.7024, 3.7223, 3.8078, 3.808, 3.8125, 3.8339, 3.8605, 3.8658, 3.8713, 3.8919, 3.9118, 3.9612, 3.9998, 4.0246, 4.0548, 4.0558, 4.1904, 4.2288, 4.3213, 4.3842, 4.4425, 4.5236, 4.5913, 4.6208, 4.6848, 4.6911, 4.6916, 4.7735, 4.7874, 4.794, 4.8415, 4.8639, 4.871, 4.9216, 4.9305, 4.932, 4.9491, 5.0137, 5.0658, 5.066, 5.0807, 5.271, 5.3019, 5.3891, 5.3942, 5.4995, 5.7012, 5.7433, 5.7482, 5.7803, 5.7886, 5.8666, 5.868, 6.025, 6.1046, 6.1366, 6.2138, 6.2556, 6.2578, 6.2948, 6.3285, 6.4076, 6.4339, 6.6006, 6.6189, 6.6641, 6.7867, 6.822, 6.8455, 7.0006, 7.1658, 7.2079, 7.208, 7.2602, 7.2829, 7.5101, 7.5113, 7.57, 7.6057, 7.6087, 7.7833, 7.9494, 7.9787, 8.0451, 8.0859, 8.1161, 8.1361, 8.1395, 8.3796, 8.4503, 8.5977, 8.6304, 8.6412, 8.771, 8.9343, 8.9413, 8.979, 9.0034, 9.0179, 9.1283, 9.2646, 9.2819, 9.46, 9.6075, 10.0297, 10.4142, 10.4334, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11), 
           "A" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1), 
           "X" = c(0.591, -0.7045, 0.1998, -1.2, -0.2748, 0.6488, -0.0464, 0.0709, -1.2, -0.1084, 0.1463, 1.2, 0.7424, 1.2, 1.1246, -0.2878, 0.6152, -0.42, 1.2, 1.2, 0.3618, 1.1865, 0.0818, -1.2, 0.8593, -0.1021, 1.1739, 0.2807, -0.8726, 1.0116, 0.4228, 0.6142, -0.3044, -0.7147, -0.7772, -0.5083, -0.3825, 0.7667, 1.2, 0.5253, 1.2, 1.2, -0.4935, -1.2, 0.7212, -0.9438, 0.2387, -1.2, 0.5291, -0.1674, 1.1412, -0.9328, 0.6634, 0.3333, -0.7123, 1.2, -0.4832, 0.3565, 0.5722, -1.0329, -0.2241, -0.8411, 1.2, 0.7209, 0.6254, 0.1107, -0.2716, 0.0051, -1.0889, 1.2, 1.2, -0.423, 1.2, 0.4794, -0.4668, 0.2772, -0.5024, -0.3908, -0.1373, 0.6094, 0.1762, 0.8822, 0.1645, -0.6743, -0.5285, -0.3291, 1.0172, 0.0233, -0.4747, 1.2, -1.2, 0.0924, -1.2, 0.3483, -0.386, 0.5705, 1.2, 1.1889, -0.1459, 1.1521, 1.2, 0.2255, 0.3866, 0.1753, -1.2, 1.2, -0.6759, 0.6161, -0.0756, -0.1158, -0.0654, 1.2, -0.0816, -0.8717, 0.3017, -0.982, -0.4491, 0.1316, 1.05, -0.9845, 0.3785, -0.8026, 0.2651, -0.9233, 1.2, -0.1918, 1.2, -0.7468, 0.7581, -0.3449, -1.2, -0.079, 0.0295, -0.4644, 0.1143, 0.8306, 0.0565, 1.2, -0.5819, 1.2, 0.572, -0.5393, -0.1491, 0.4831, 0.0408, -0.9452, 0.6162, 0.8923, 1.2, -0.0585, 0.0334, -0.6505, -1.0751, 0.6102, 0.717, -0.6497, -1.2, -0.3648, 0.4697, -1.0421, 0.1834, -1.2, -0.7025, -0.1347, -0.785, -0.8318, -0.7814, -0.2372, -0.2614, -1.0397, -1.2, 0.6774, -0.8028, -1.0292, -1.2, 0.0708, -0.0817, -0.4429, 0.413, -0.1401, -0.8748, -1.2, -0.7557, -1.0771, -0.838, -0.446, -0.8366, -0.526, -0.5437, -1.2, -0.6688, 0.4963, -1.2, -0.495, 0.779, -1.2, -0.7388, -0.6316, -0.6253, -0.564, 0.5762, 1.1168, 0.0583, 0.2203, -0.1059, 0.6471, -0.0811, 0.6312, 0.0156, -0.8045, -1.1241, 0.4943, -0.2819, -0.1622, 0.1603, -1.2, 0.5105, 0.4579, -0.7416, 0.8165, -1.193, -0.5748, -0.424, -1.0011, -1.2, -1.0157, 0.688, -0.8993, -0.2307, -0.0179, -0.1043, -0.8387, 0.2837, 0.1785, -0.6486, -1.2, -0.399, -0.8616, -1.2, -0.4967, -1.2, -0.7369, -5e-04, 0.1076, -1.0822, -0.2976, 0.617, -1.2, -0.428, -1.2))
d$A2 <- as.factor(d$A)
d$X2 <- d$X^2
d$Y1 <- as.numeric(d$time<=5)*as.numeric(d$status==1)

## *** prepare
## create weights for AIPTW
e.T <- glm(A2~X,family="binomial",data=d)
myW.T1 <- (d$A2=="1")/predict(e.T, type = "response")
myW.T0 <- (d$A2=="0")/(1-predict(e.T, type = "response"))
## range(1/d$PS2 - myW.T)
## range((iW.IPTW[,2]-iW.IPTW[,1]) - myW.T)
e.C <- coxph(Surv(time,status==0)~1,data=d, x = TRUE)
myW.C <- (1-(d$time < 5)*(d$status==0))/predictCoxPL(e.C, newdata = d, times = pmin(d$time,5), diag = TRUE)$survival
## range(iW.IPCW - myW.C)
## range(d[["IPCW"]] - myW.C)

## create counterfactual data
d1 <- copy(d)
d1$A <- 1
d1$A2 <- factor(1,levels = 0:1)
d0 <- copy(d)
d0$A <- 0
d0$A2 <- factor(0,levels = 0:1)

## *** test
test_that("[ate] IPCW LR - censoring (based on Paul's script and results)", {
    
    ## no censoring: G-formula on glm 
    e.glm <- glm(Ytrue~X+A2,family="binomial",data=d)
    eATE.wglm <- ate(e.glm, data = d, treatment = "A2", band = FALSE, verbose = FALSE)

    expect_equal(mean(predictRisk(e.glm, newdata = d0)),
                 0.2468305, tol = 1e-5)
    expect_equal(mean(predictRisk(e.glm, newdata = d1)),
                 0.4593332, tol = 1e-5)
    expect_equal(eATE.wglm$meanRisk$meanRisk.Gformula,
                 c(0.2468305,0.4593332), tol = 1e-5)

    ## censoring:
    ## IPCW model
    suppressWarnings(eW.glm <- glm(Y1~X2+A2,family="binomial",data=d, weights = myW.C))

    e.wglm <- wglm(regressor.event = ~ X2+A2,
                   formula.censor = Surv(time,status==0) ~ 1,
                   times = 5, data = d, product.limit = TRUE)

    
    ## G-formula on IPCW glm
    eATE.wglm <- ate(e.wglm, data = d, treatment = "A2", times = 5, band = FALSE, verbose = FALSE)

    expect_equal(mean(predictRisk(eW.glm, newdata = d0)), 0.1886017, tol = 1e-5) 
    expect_equal(mean(predictRisk(eW.glm, newdata = d1)), 0.4925787, tol = 1e-5)
    ## minor difference due to coxph using Newton Raphson instead of explicit formulae
    
    expect_equal(eATE.wglm$meanRisk$meanRisk.Gformula,
                 c(0.1886017,0.4925787), tol = 1e-3)
    ## minor difference due to difference weights computations

    ## AIPTW on IPCW glm [:TOFIX:]
    ## eATE.wglm2 <- ate(e.wglm, treatment = A2 ~ X, censor = Surv(time,status==0) ~ 1,
    ##                   product.limit = TRUE, estimator = c("Gformula","AIPTW"),
    ##                   data = d, times = 5, band = FALSE, verbose = FALSE)

    ## expect_equal(eATE.wglm2$meanRisk[1,meanRisk.AIPTW],
    ##              mean( myW.T0*myW.C*(d$Y1 - fitted(eW.glm)) + predict(eW.glm, newdata = d0, type = "response") ),
    ##              tol = 1e-4)
                 
    ## expect_equal(eATE.wglm2$meanRisk[2,meanRisk.AIPTW],
    ##              mean( myW.T1*myW.C*(d$Y1 - fitted(eW.glm)) + predict(eW.glm, newdata = d1, type = "response") ),
    ##              tol = 1e-4)
                 
    ## mean( myW.T1*myW.C*(d$Y1 - fitted(eW.glm)) + predict(eW.glm, newdata = d1, type = "response") )

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
                  band = FALSE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate1 <- as.data.table(confint(e.ate1, band = TRUE, p.value = TRUE))
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose, cause = 1,
                  iid = TRUE, product.limit = FALSE
                  )
    set.seed(15)
    dt.ate2 <- as.data.table(confint(e.ate2, band = TRUE, p.value = TRUE))

    expect_equal(dt.ate1[,.SD,.SDcols = names(dt.ate2)],  dt.ate2)

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
        iDT[, iid.Gformula := iid.Gformula + rowMeans(iPred.risk$absRisk.iid[,1,])]

        iDT[, iid.IPTW := (X1test*Y/pi - mean(X1test*Y/pi))/.N]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y/pi^2))]

        iDT[, iid.AIPTW := (X1test*Y/pi + r*(1-X1test/pi) - mean(X1test*Y/pi + r*(1-X1test/pi)))/.N]
        iDT[, iid.AIPTW := iid.AIPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*(Y-r)/pi^2))]
        iDT[, iid.AIPTW := iid.AIPTW + rowMeans(rowMultiply_cpp(iPred.risk$absRisk.iid[,1,], scale = 1-X1test/pi))]

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
    GS <- data.table("value" = c(0.2794404, 0.1838809, 0.3965070), 
                     "se" = c(0.2272017, 0.1549114, 0.2688568), 
                     "lower" = c(-0.1658668, -0.1197399, -0.1304426), 
                     "upper" = c(0.7247475, 0.4875017, 0.9234565), 
                     "p.value" = c(0.2187263, 0.2352249, 0.1402693), 
                     "quantileBand" = c(1.969141, 1.922815, 1.933780), 
                     "adj.p.value" = c(0.2235, 0.2291, 0.1376), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))

    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(3.999770, 3.187715, 6.052789), 
                     "se" = c(2.266232, 2.180114, 3.865633), 
                     "lower" = c(0, 0, 0), 
                     "upper" = c( 8.441504,  7.460660, 13.629289), 
                     "p.value" = c(0.1856088, 0.3156261, 0.1911770), 
                     "quantileBand" = c(1.971403, 1.958227, 1.974396), 
                     "adj.p.value" = c(0.1872, 0.3072, 0.1973), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
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
                  band = FALSE, product.limit = FALSE, cause = 1
                  )
    set.seed(15)
    dt.ate1 <- as.data.table(confint(e.ate1, band = TRUE, p.value = TRUE))
    e.ate2 <- ate(data = dtS, times = tau,
                  event = e.S,
                  censor = e.C,
                  treatment = e.T,
                  method.iid = 2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  band = FALSE, product.limit = FALSE, cause = 1
                  )
    set.seed(15)
    dt.ate2 <- as.data.table(confint(e.ate2, band = TRUE, p.value = TRUE))

    expect_equal(dt.ate1, dt.ate2, tol = 1e-8)
    
    ## ## manual calculation ## ##
    iPred.Cens <- predictCox(e.C, newdata = dtS, times = pmin(tau,dtS$time-1e-10), type = "survival", iid = TRUE)
    iPred.Cens$survivalDiag <- diag(iPred.Cens$survival)
    iPred.Cens$survivalDiag.iid <- sapply(1:n, function(iObs){iPred.Cens$survival.iid[,iObs,iObs]})

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
        iDT[, iid.Gformula := iid.Gformula + rowMeans(iPred.risk_tau$absRisk.iid[,1,])]

        ## IPTW
        iDT[, iid.IPTW := (X1test*Y*C/(pi*G) - mean(X1test*Y*C/(pi*G)))/.N]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -X1test*Y*C/(pi^2*G)))]
        iDT[, iid.IPTW := iid.IPTW + rowMeans(rowMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]

        ## AIPTW,AIPCW ##        
        iDT[, iid.AIPTW.ate := (X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi - mean(X1test*Y*C/(pi*G) + r*(1-X1test/pi) + I*X1test/pi))/.N]

        ## outcome
        iDT[, iid.AIPTW.outcome := rowMeans(rowMultiply_cpp(iPred.risk_tau$absRisk.iid[,1,], (1-X1test/pi)))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + rowMeans(rowMultiply_cpp(augTermY1.iid,X1test/pi))]
        iDT[, iid.AIPTW.outcome := iid.AIPTW.outcome + rowMeans(rowMultiply_cpp(augTermY2.iid,X1test/pi))]

        ## treatment
        iDT[, iid.AIPTW.treatment := rowMeans(rowMultiply_cpp(attr(iPred.logit,"iid"), scale = -(X1test/pi^2)*(Y*C/G - r + I)))]

        ## survival
        iDT[, iid.AIPTW.survival := rowMeans(rowMultiply_cpp(augTermSurv.iid, X1test/pi))]

        ## censoring
        iDT[, iid.AIPTW.censoring := rowMeans(rowMultiply_cpp(iPred.Cens$survivalDiag.iid, scale = -X1test*Y*C/(pi*G^2)))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + rowMeans(rowMultiply_cpp(augTermG1.iid,X1test/pi))]
        iDT[, iid.AIPTW.censoring := iid.AIPTW.censoring + rowMeans(rowMultiply_cpp(augTermG2.iid,X1test/pi))]

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
                     "quantileBand" = c(1.958676, 1.963257, 1.952228), 
                     "adj.p.value" = c(0.1177, 0.2880, 0.0331), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "diffAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
                 tol = 1e-6)

    GS <- data.table("value" = c(3.102436, 1.555969, 3.371188), 
                     "se" = c(1.4694536, 0.5500765, 1.2732042), 
                     "lower" = c(0.2223595, 0.4778388, 0.8757540), 
                     "upper" = c(5.982512, 2.634099, 5.866623), 
                     "p.value" = c(0.15249898, 0.31215422, 0.06254973), 
                     "quantileBand" = c(1.955519, 1.941080, 1.992390), 
                     "adj.p.value" = c(0.1553, 0.3178, 0.0605), 
                     "estimator" = c("Gformula", "IPTW", "AIPTW"))
    expect_equal(GS, dt.ate1[type == "ratioAte",.(value,se,lower,upper,p.value,quantileBand,adj.p.value,estimator)],
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
## * [ate] Case only
cat("[ate] Case only \n")
set.seed(11)
d <- sampleData(100, outcome = "survival")
d$id <- 1:NROW(d)
tau <- 3

## ** statistical model
e.cox <- coxph(Surv(time,event) ~ X1 + X2, x = TRUE, y = TRUE, data = d)

## ** ate: case only
test_that("Case only ate", {

    ## automatically
    e.ateC <- ate(e.cox,
                  treatment  = "X1",
                  data = d[X1 == "1"], time = tau, se = TRUE,
                  data.index = d[X1 == "1",id],
                  verbose = 0)

    ## manually
    d0C <- copy(d[X1 == "1"])
    d0C[, X1:= factor("0", levels = c("0","1"))]
    d1C <- copy(d[X1 == "1"])
    d1C[, X1:= factor("1", levels = c("0","1"))]

    e.pred0C <- predictCox(e.cox, newdata = d0C, times = tau, average.iid = TRUE)
    e.pred1C <- predictCox(e.cox, newdata = d1C, times = tau, average.iid = TRUE)
    e2.ateC <- c(mean(1-e.pred0C$survival),
                 mean(1-e.pred1C$survival),
                 mean( (1-e.pred1C$survival)-(1-e.pred0C$survival)),
                 mean(1-e.pred1C$survival)/mean(1-e.pred0C$survival))

    iid0C <- -e.pred0C$survival.average.iid
    iid0C[d0C$id,] <- iid0C[d0C$id,] + (1-e.pred0C$survival-e2.ateC[1])/NROW(d0C)
    iid1C <- -e.pred1C$survival.average.iid
    iid1C[d1C$id,] <- iid1C[d1C$id,] + (1-e.pred1C$survival-e2.ateC[2])/NROW(d1C)

    e2.se.ateC <- c(sqrt(sum(iid0C^2)),
                    sqrt(sum(iid1C^2)),
                    sqrt(sum((iid1C-iid0C)^2)),
                    sqrt(sum((iid1C/e2.ateC[1]-e2.ateC[2]*iid0C/e2.ateC[1]^2)^2))
                    )

    
    expect_equal(as.data.table(e.ateC)$value, e2.ateC, tol = 1e-6)
    expect_equal(as.data.table(e.ateC)$se, e2.se.ateC, tol = 1e-5)
})
## * [ate] Bootstrap
cat("[ate] Bootstrap \n")
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

    e.ate1 <- ate(data = dtS, times = tau,
                  event = e.S,
                  treatment = e.T,
                  censor = e.C,
                  verbose = verbose,
                  se = TRUE, iid = FALSE, product.limit = FALSE,
                  B = 5, cpus = 1
                  )
    xx <- capture.output(print(e.ate1))

    e.ate2 <- ate(data = dtS, times = tau,
                  event = Surv(time, event) ~ X1f + strat(X2) + X3*X6,
                  treatment = X1f ~ X2,
                  censor = Surv(time, event == 0) ~ X1f + X2,
                  estimator = c("Gformula","IPTW","AIPTW"),
                  verbose = verbose,
                  se = TRUE, iid = FALSE, product.limit = FALSE,
                  B = 5, cpus = 1
                  )
    yy <- capture.output(print(e.ate2))

})

## * [ate] Miscellaneous
cat("[ate] Miscellaneous \n")

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
## ** Ate using list of formulae
test_that("[ate] using list of formulae",{
    set.seed(10)
    d <- sampleData(100)
    GS <- ate(list(Hist(time,event)~X1+X2,
                   Hist(time,event)~X1+X2),
              treatment = "X1", data = d, times = 1:5, cause = 1,
              verbose = verbose)
})

## ** Specification of the censoring mecanism
test_that("[ate] Specification of the censoring mecanism",{
    set.seed(10)
    d <- sampleData(100)
    d$allevent <- as.numeric(d$event>0)
        
    expect_error(ate(event = Surv(time,allevent)~X1+X2,
                     censor = Surv(time,allevent)~X1,
                     treatment = X1 ~ X2,
                     data = d, times = 1:5, cause = 1,
                     verbose = verbose))

    expect_error(ate(event = Hist(time,event)~X1+X2,
                     censor = Surv(time,allevent)~X1,
                     treatment = X1 ~ X2,
                     data = d, times = 1:5, cause = 1,
                     verbose = verbose))

    ## should run
    test1 <- ate(event = Hist(time,event)~X1+X2,
                 censor = Surv(time,allevent==0)~X1,
                 treatment = X1 ~ X2,
                 data = d, times = 1:5, cause = 1,
                 verbose = verbose)
})


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

## ** Before the first jump or at the first jump
test_that("[ate] double robust estimator - before or at the first jump",{
    ## Error in rowSums(integrand.St[[iGrid]][, 1:beforeTau.nJumpC[iTau]]) : 
    ##   'x' must be an array of at least two dimensions

    set.seed(1)
    d <- sampleData(100,outcome="survival")
    tau <- c(d[event>0,min(time)] + c(-1e-5,+1e-5), median(d$time))
    
    fit <-  ate(event = coxph(Surv(time,event)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE),
                treatment = X1~X4,
                censor = coxph(Surv(time,event==0)~X1+X4+X7+X8,data=d,x=TRUE,y=TRUE),
                data=d, time=tau, se = TRUE, verbose = verbose)
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

    GS <- data.table("type" = c("ate", "ate", "ate", "ate", "diffAte", "diffAte", "ratioAte", "ratioAte"), 
                     "level" = c("0", "1", "0", "1", "0.1", "0.1", "0.1", "0.1"), 
                     "time" = c(3, 3, 4, 4, 3, 4, 3, 4), 
                     "value" = c(0.2585028, 0.5034093, 0.3052605, 0.5005653, 0.2449065, 0.1953049, 1.9474034, 1.6397974), 
                     "se" = c(0.04636428, 0.15684698, 0.04902926, 0.15630839, 0.16327935, 0.16374745, 0.69821732, 0.57546242), 
                     "lower" = c( 0.16763051,  0.19599485,  0.20916489,  0.19420653, -0.07511518, -0.12563425,  0.57892260,  0.51191176), 
                     "upper" = c(0.3493751, 0.8108237, 0.4013561, 0.8069241, 0.5649281, 0.5162440, 3.3158842, 2.7676830), 
                     "p.value" = c(NA, NA, NA, NA, 0.1336343, 0.2329791, 0.1748165, 0.2662254), 
                     "estimator" = c("AIPTW", "AIPTW", "AIPTW", "AIPTW", "AIPTW", "AIPTW", "AIPTW", "AIPTW"))
    
    expect_equal(as.data.table(e.ateRR), GS, tol = 1e-6)
    
    
    e.ateRR <- ate(event = CSC(Hist(time,event) ~ X1 + X2 + X3, data = dtS),
                   treatment = "X1",
                   data = dtS, times = 3:4, verbose = verbose, cause = 1
                   )

})

## ** Tied events
## Brice: Tuesday, Sept 1, 2020 11:10
n <- 500
set.seed(10)
dt <- sampleData(n, outcome="survival")
dt$time <- round(dt$time,1)+0.1
dt$X1 <- factor(rbinom(n, prob = c(0.4) , size = 1), labels = paste0("T",0:1))

test_that("ate in presence of ties",{
    m.event <-  coxph(Surv(time,event)~ X1,data=dt,x=TRUE,y=TRUE)

    ateRobust <- ate(event = m.event,
                     treatment = "X1",
                     data = dt, times = 1, 
                     cause = 1, product.limit = FALSE,
                     verbose = verbose)

    iPred <- predictCox(m.event, newdata = unique(dt[,.SD,.SDcols = "X1"]),
                        times = 1, se = TRUE, type = "survival", keep.newdata = TRUE)

    iDT <- as.data.table(iPred)
    setkeyv(iDT, "X1")
    expect_equal(ateRobust$meanRisk$meanRisk.Gformula, 1-iDT$survival, tol = 1e-7)
    expect_equal(ateRobust$meanRisk$meanRisk.Gformula.se, iDT$survival.se, tol = 1e-7)
})
