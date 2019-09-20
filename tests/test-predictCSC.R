### test-predictCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: sep 18 2019 (14:18) 
##           By: Brice Ozenne
##     Update #: 244
#----------------------------------------------------------------------
## 
### Commentary: 
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Settings
library(riskRegression)
library(testthat)
library(mstate)
tmat <- trans.comprisk(2, names = c("0", "1", "2"))
library(survival)
library(timereg); nsim.band <- 100

## * [predictCSC] predict.CSC vs. manual calculation and mstate on small examples
cat("[predictCSC] predict.CSC vs. manual calculation on small examples \n")
## ** Data
## simulatenous events
df1 <- data.frame(time = rep(1:10,2),
                 event = c(rep(1,10),rep(2,10)),
                 X1 = 0:1
                 )
df1$event1 <- as.numeric(df1$event == 1)
df1$event2 <- as.numeric(df1$event == 2)

set.seed(11)
dfS <- rbind(cbind(df1, grp = 1, X2 = rbinom(20,1,.4)),
             cbind(rbind(df1,df1), grp = 2, X2 = rbinom(20,1,.4)),
             cbind(df1, grp = 3, X2 = rbinom(20,1,.4))
             )

## distinct events
df2 <- data.frame(time = c(1:10,-0.1+1:10),
                  event = c(rep(1,10),rep(2,10)),
                  X1 = 0:1
                  )
df2$event1 <- as.numeric(df2$event == 1)
df2$event2 <- as.numeric(df2$event == 2)

## ** No risk factor no strata
test_that("[predictCSC] (no strata, data=df1): compare to manual estimation",{

    CSC.exp <- CSC(Hist(time,event) ~ 1, data = df1)
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = FALSE)  
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE) 
    CSC.prodlim <- CSC(Hist(time,event) ~ 1, data = df1, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    expect_equal(pred.prodlim,pred.exp2)
    CSC.prodlimE <- CSC(Hist(time,event) ~ 1, data = df1, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  

    ## test baseline hazard
    lambda1 <- lambda2 <- 1/seq(20,2,by = -2)
    lambda2E <- predictCox(CSC.prodlimE$models[[2]], type = "hazard")$hazard
    expect_equal(predictCox(CSC.exp$models[[1]], type = "hazard")$hazard,
                 lambda1)
    expect_equal(predictCox(CSC.exp$models[[2]], type = "hazard")$hazard,
                 lambda2)
    expect_equal(predictCox(CSC.prodlim$models[[1]], type = "hazard")$hazard,
                 lambda1)
    expect_equal(predictCox(CSC.prodlim$models[[2]], type = "hazard")$hazard,
                 lambda1+lambda2)
    expect_equal(basehaz(CSC.prodlimE$models[[2]])$hazard,
                 cumsum(lambda2E)
                 )

    ## test absolute risk
    survival <- c(1,exp(cumsum(-lambda1-lambda2)))[1:10]
    expect_equal(as.double(pred.exp$absRisk),
                 cumsum(lambda1*survival))
    survival <- c(1,cumprod(1-(lambda1+lambda2)))[1:10]
    expect_equal(as.double(pred.prodlim$absRisk),
                 cumsum(lambda1*survival))
    survival <- c(1,cumprod(1-lambda2E))[1:10]
    expect_equal(as.double(pred.prodlimE$absRisk),
                 cumsum(lambda1*survival))
})

test_that("[predictCSC] (no strata, data=df1): compare to manual estimation",{
    CSC.exp <- CSC(Hist(time,event) ~ 1, data = df2)
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = FALSE)  
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    CSC.prodlim <- CSC(Hist(time,event) ~ 1, data = df2, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    expect_equal(pred.exp2, pred.prodlim)
    CSC.prodlimE <- CSC(Hist(time,event) ~ 1, data = df2, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    ## test baseline hazard
    lambda1 <- 1/seq(19,1,by = -2)
    lambda2 <- 1/seq(20,2,by = -2)
    lambda2E <- as.data.table(predictCox(CSC.prodlimE$models[[2]], type = "hazard"))
    expect_equal(predictCox(CSC.exp$models[[1]], times = 1:10, type = "hazard")$hazard,
                 lambda1)
    expect_equal(predictCox(CSC.exp$models[[2]], times = -0.1 + 1:10, type = "hazard")$hazard,
                 lambda2)
    expect_equal(predictCox(CSC.prodlim$models[[1]], times = 1:10, type = "hazard")$hazard,
                 lambda1)
    expect_equal(predictCox(CSC.prodlim$models[[2]], type = "hazard")$hazard,
                 as.double(rbind(lambda2,lambda1)))
    expect_equal(basehaz(CSC.prodlimE$models[[2]])$hazard,
                 cumsum(lambda2E$hazard)
                 )
    ## test absolute risk
    vec.lambda <- cumsum(as.double(rbind(lambda2,lambda1)))
    survival <- exp(-matrix(vec.lambda, ncol = 2, byrow = TRUE)[,1])
    expect_equal(as.double(pred.exp$absRisk),
                 cumsum(lambda1*survival))
    vec.lambda <- as.double(rbind(lambda2,lambda1))
    survival <- matrix(cumprod(1-vec.lambda), ncol = 2, byrow = TRUE)[,1]
    expect_equal(as.double(pred.prodlim$absRisk),
                 cumsum(lambda1*survival))
    vec.lambda <-lambda2E$hazard
    survival <- matrix(cumprod(1-vec.lambda), ncol = 2, byrow = TRUE)[,1]
    expect_equal(as.double(pred.prodlimE$absRisk),
                 cumsum(lambda1*survival))
})

test_that("[predictCSC]: compare to mstate",{
    for(iData in 1:2){# iData <- 1
        df <- switch(as.character(iData),
                     "1"=df1,
                     "2"=df2)
        ## fit CSC
        CSC.exp <- CSC(Hist(time,event) ~ 1, data = df)
        ## fit mstate
        dL <- msprep(time = c(NA, "time", "time"),
                     status = c(NA,"event1", "event2"),
                     data = df, keep = c("X1"),
                     trans = tmat)
        dL.exp <- expand.covs(dL,  c("X1"))
        e.coxph <- coxph(Surv(time, status) ~ strata(trans),
                         data = dL.exp)
        ## compare
        newdata.L <- data.frame(trans = c(1, 2), strata = c(1, 2))
        newdata <- data.frame(NA)
        pred.msfit <- msfit(e.coxph, newdata = newdata.L, trans = tmat)
        suppressWarnings(
            pred.probtrans <- probtrans(pred.msfit,0)[[1]]
        )
        pred.RR <- predict(CSC.exp,
                           times = pred.probtrans[,"time"],
                           newdata = newdata,
                           cause = 1,
                           se = FALSE,
                           product.limit = TRUE)
        expect_equal(pred.probtrans[,"pstate2"],
                     as.double(pred.RR$absRisk)
                     )
    }
})


test_that("[predict.CSC] check absolute risks add up to one",{
    for(iData in 1:2){
        df <- switch(as.character(iData),
                     "1"=df1,
                     "2"=df2)
        ## fit CSC
        CSC.exp <- CSC(Hist(time,event) ~ 1, data = df, method = "breslow")
        ## compute all transition probabilites
        seqTimes <- sort(unique(df$time))
        nTimes <- length(seqTimes)
        newdata <- data.frame(NA)
        hazAll <- predictCox(CSC.exp$models[[1]], times = seqTimes, newdata = newdata, type = "hazard")$hazard+predictCox(CSC.exp$models[[2]], times = seqTimes, newdata = newdata, type = "hazard")$hazard
        surv.prodlim <- cumprod(1-as.double(hazAll))
        haz1.prodlim <- predict(CSC.exp, times = seqTimes, newdata = newdata, cause = 1, product.limit = TRUE)  
        haz2.prodlim <- predict(CSC.exp, times = seqTimes, newdata = newdata, cause = 2, product.limit = TRUE)  
        expect_equal(surv.prodlim+as.double(haz1.prodlim$absRisk)+as.double(haz2.prodlim$absRisk),rep(1,nTimes))
    }
})

## ** No risk factor strata
test_that("[predict.CSC] (strata): compare to manual estimation",{
    ## as.data.table(dfS)[,.N, by = "grp"]
    CSC.exp <- CSC(Hist(time,event) ~ strata(grp), data = dfS, method = "breslow")    
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(grp = 1:3), cause = 1, product.limit = FALSE)
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(grp = 1:3), cause = 1, product.limit = TRUE)  

    CSC.prodlim <- CSC(Hist(time,event) ~ strata(grp), data = dfS, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(grp = 1:3), cause = 1, product.limit = TRUE)  
    expect_equal(pred.exp2,pred.prodlim)
    
    CSC.prodlimE <- CSC(Hist(time,event) ~ strata(grp), data = dfS, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(grp = 1:3), cause = 1, product.limit = TRUE)  

    ## baseline
    baseline1 <- as.data.table(predictCox(CSC.exp$models[[1]], type = "hazard"))
    lambda1 <- baseline1[strata=="grp=1",hazard]
    expect_equal(lambda1, 1/seq(20,2,by = -2))

    baseline2 <- as.data.table(predictCox(CSC.exp$models[[1]], type = "hazard"))
    lambda2 <- baseline2[strata=="grp=1",hazard]
    expect_equal(lambda2, 1/seq(20,2,by = -2))
        
    baseline1E <- as.data.table(predictCox(CSC.prodlimE$models[[1]], type = "hazard"))
    lambda1E.1 <- baseline1E[strata=="grp=1",hazard]
    lambda1E.2 <- baseline1E[strata=="grp=2",hazard]
    expect_equal(lambda1E.1, 1/seq(20,2,by = -2))

    baseline2E <- as.data.table(predictCox(CSC.prodlimE$models[[2]], type = "hazard"))
    lambda2E.1 <- baseline2E[strata=="grp=1",hazard]
    lambda2E.2 <- baseline2E[strata=="grp=2",hazard]
    ## 
        
    ## test absolute risk
    survival <- c(1,exp(cumsum(-lambda1-lambda2)))[1:10]
    expect_equal(as.double(pred.exp$absRisk[1,]),
                 cumsum(lambda1*survival))
    expect_equal(as.double(pred.exp$absRisk[1,]),
                 as.double(pred.exp$absRisk[2,]))
    expect_equal(as.double(pred.exp$absRisk[1,]),
                 as.double(pred.exp$absRisk[3,]))

    survival <- c(1,cumprod(1-(lambda1+lambda2)))[1:10]
    expect_equal(as.double(pred.prodlim$absRisk[1,]),
                 cumsum(lambda1*survival))
    expect_equal(as.double(pred.prodlim$absRisk[1,]),
                 as.double(pred.prodlim$absRisk[2,]))
    expect_equal(as.double(pred.prodlim$absRisk[1,]),
                 as.double(pred.prodlim$absRisk[3,]))

    survival <- c(1,cumprod(1-lambda2E.1))[1:10]
    expect_equal(as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=1"]),
                 cumsum(lambda1E.1*survival))
    survival <- c(1,cumprod(1-lambda2E.2))[1:10]
    expect_equal(as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=2"]),
                 cumsum(lambda1E.2*survival))    
    expect_equal(as.double(pred.prodlimE$absRisk[1,]),
                 as.double(pred.prodlimE$absRisk[3,]))
})

## ** risk factor no strata
test_that("[predict.CSC] (risk factor, strata): compare to manual estimation",{
    CSC.exp <- CSC(Hist(time,event) ~ X1, data = df1)
    CSC.prodlim <- CSC(Hist(time,event) ~ X1, data = df1, surv.type = "survival", method = "breslow")
    CSC.prodlimE <- CSC(Hist(time,event) ~ X1, data = df1, surv.type = "survival", method = "efron")

    calcCIF <- function(CSC, X){
        lambda01 <- as.data.table(predictCox(CSC$models[[1]], centered = FALSE, type = "hazard"))
        lambda02 <- as.data.table(predictCox(CSC$models[[2]], centered = FALSE, type = "hazard"))
        eXb.1 <- as.double(exp(X%*%coef(CSC)[[1]]))
        eXb.2 <- as.double(exp(X%*%coef(CSC)[[2]]))
        if(CSC$surv.type == "survival"){
            survival <- cumprod(1-eXb.2*lambda02$hazard)
        }else{           
            Lambda01 <- cumsum(lambda01$hazard)
            Lambda02 <- cumsum(lambda02$hazard)
            survival <- exp(-eXb.1*Lambda01-eXb.2*Lambda02)
        }
        survival <- c(1,survival[-length(survival)])
        absRisk <- cumsum(eXb.1*lambda01$hazard*survival)
        return(absRisk)
    }
    for(iX1 in 0:1){
        newdata <- data.frame(X1=iX1)
        pred.exp <- predict(CSC.exp, times = 1:10, newdata = newdata, cause = 1, product.limit = FALSE)  
        expect_equal(as.double(pred.exp$absRisk),
                     calcCIF(CSC.exp, X = as.matrix(newdata))
                     )
        pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = newdata, cause = 1, product.limit = TRUE)
        expect_equal(as.double(pred.prodlim$absRisk),
                     calcCIF(CSC.prodlim, X = as.matrix(newdata))
                     )
        pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(as.double(pred.prodlimE$absRisk),
                     calcCIF(CSC.prodlimE, X = as.matrix(newdata))
                     )
    }
})

test_that("[predict.CSC] - compare to manual estimation",{
    CSC.exp <- CSC(Hist(time,event) ~ X1, data = df2, method = "breslow")

    CSC.prodlim <- CSC(Hist(time,event) ~ X1, data = df2, surv.type = "survival", method = "breslow")
    
    CSC.prodlimE <- CSC(Hist(time,event) ~ X1, data = df2, surv.type = "survival", method = "efron")

    calcCIF <- function(CSC, X){

        lambda01 <- as.data.table(predictCox(CSC$models[[1]], centered = FALSE, type = "hazard"))
        lambda02 <- as.data.table(predictCox(CSC$models[[2]], centered = FALSE, type = "hazard"))
        eXb.1 <- as.double(exp(X%*%coef(CSC)[[1]]))
        eXb.2 <- as.double(exp(X%*%coef(CSC)[[2]]))

        if(CSC$surv.type == "survival"){
            survival <- cumprod(1-eXb.2*lambda02$hazard)
        }else{           
            Lambda01 <- cumsum(lambda01$hazard)
            Lambda02 <- cumsum(lambda02$hazard)
            survival <- exp(-eXb.1*Lambda01-eXb.2*Lambda02)
        }

        survival <- c(1,survival[-length(survival)])
        absRisk <- cumsum(eXb.1*lambda01$hazard*survival)
        return(absRisk)
        
    }
 
    for(iX1 in 0:1){ # iX1 <- 0
        newdata <- data.frame(X1=iX1)
        pred.exp <- predict(CSC.exp, times = sort(unique(df2$time)), newdata = newdata, cause = 1, product.limit = FALSE)  
        expect_equal(as.double(pred.exp$absRisk),
                     calcCIF(CSC.exp, X = as.matrix(newdata))
                     )
    
        pred.prodlim <- predict(CSC.prodlim, times = sort(unique(df2$time)), newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(as.double(pred.prodlim$absRisk),
                     calcCIF(CSC.prodlim, X = as.matrix(newdata))
                     )

        pred.prodlimE <- predict(CSC.prodlimE, times = sort(unique(df2$time)), newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(as.double(pred.prodlimE$absRisk),
                     calcCIF(CSC.prodlimE, X = as.matrix(newdata))
                     )
    }
})


test_that("[predict.CSC]: compare to mstate",{
    for(iData in 1:2){
        df <- switch(as.character(iData),
                     "1"=df1,
                     "2"=df2)
        ## fit CSC
        CSC.exp <- CSC(Hist(time,event) ~ X1, data = df, method = "breslow")
        ## fit mstate
        df$event1 <- as.numeric(df$event == 1)
        df$event2 <- as.numeric(df$event == 2)
        dL <- msprep(time = c(NA, "time", "time"),
                     status = c(NA,"event1", "event2"),
                     data = df, keep = c("X1"),
                     trans = tmat)
        dL.exp <- expand.covs(dL,  c("X1"))
        e.coxph <- coxph(Surv(time, status) ~ X1.1 + X1.2 + strata(trans),
                         data = dL.exp)
        ## compare
        for(iX1 in 0:1){ # iX1 <- 0
            newdata.L <- data.frame(X1.1 = c(iX1,0), X1.2 = c(0,iX1), trans = c(1, 2), strata = c(1, 2))    
            newdata <- data.frame(X1 = iX1)
            pred.msfit <- msfit(e.coxph, newdata = newdata.L, trans = tmat)
            suppressWarnings(
                pred.probtrans <- probtrans(pred.msfit,0)[[1]]
            )
            pred.RR <- predict(CSC.exp, times = pred.probtrans[,"time"], newdata = newdata, cause = 1, product.limit = TRUE)
            expect_equal(pred.probtrans[,"pstate2"],
                         as.double(pred.RR$absRisk)
                         )
        }
    }
})

test_that("[predict.CSC]: check absolute risks add up to one",{
    for(iData in 1:2){
        df <- switch(as.character(iData),
                     "1"=df1,
                     "2"=df2)

        ## fit CSC
        CSC.exp <- CSC(Hist(time,event) ~ X1, data = df, method = "breslow")

        ## compute all transition probabilites
        seqTimes <- sort(unique(df$time))
        nTimes <- length(seqTimes)
  
        for(iX1 in 0:1){            
            newdata <- data.frame(X1 = iX1)

            hazAll <- predictCox(CSC.exp$models[[1]], times = seqTimes, newdata = newdata, type = "hazard")$hazard+predictCox(CSC.exp$models[[2]], times = seqTimes, newdata = newdata, type = "hazard")$hazard
            surv.prodlim <- cumprod(1-as.double(hazAll))
            haz1.prodlim <- predict(CSC.exp, times = seqTimes, newdata = newdata, cause = 1, product.limit = TRUE)  
            haz2.prodlim <- predict(CSC.exp, times = seqTimes, newdata = newdata, cause = 2, product.limit = TRUE)  

            expect_equal(surv.prodlim+as.double(haz1.prodlim$absRisk)+as.double(haz2.prodlim$absRisk),rep(1,nTimes))
        }
    }
})

## * [predictCSC] predict.CSC vs. mstate on simulated and real data
cat("[predictCSC] predict.CSC vs. mstate \n")
## ** Data
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]

d <- d[event!=0]

d[, event1 := as.numeric(event == 1)]
d[, event2 := as.numeric(event == 2)]

## ** Mstate
tmat <- trans.comprisk(2, names = c("0", "1", "2"))

dL <- msprep(time = c(NA, "time", "time"),
             status = c(NA,"event1", "event2"),
             data = d, keep = c("X1","X2","X16","Xcat2"),
             trans = tmat)
dL.exp <- expand.covs(dL,  c("X1","X2","X16","Xcat2"))

## ** No covariates
test_that("predict.CSC (no covariates): compare to mstate",{
    newdata <- data.frame(NA)
    newdata.L <- data.frame(trans = c(1, 2), strata = c(1, 2))
    # mstate
    e.coxph <- coxph(Surv(time, status) ~1 + strata(trans),
                     data = dL.exp)
    
    pred.msfit <- msfit(e.coxph, newdata = newdata.L, trans = tmat)
    pred.probtrans <- probtrans(pred.msfit,0)[[1]]

    ## riskRegression
    CSC.RR1 <- CSC(Hist(time,event)~1, data = d, method = "breslow")

    pred.RR1a <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)

    pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
    
    CSC.RR2 <- CSC(Hist(time,event)~1, data = d, surv.type = "survival", method = "breslow")
    pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                        keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                        keep.newdata = FALSE, se = TRUE, product.limit = FALSE)


    expect_equal(as.double(pred.RR1a$absRisk),pred.probtrans[,"pstate2"])
    expect_equal(as.double(pred.RR1b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)
    expect_equal(as.double(pred.RR1c$absRisk),pred.probtrans[,"pstate3"])
    expect_equal(as.double(pred.RR1d$absRisk),pred.probtrans[,"pstate3"], tol = 5e-2)
    expect_equal(as.double(pred.RR2a$absRisk),pred.probtrans[,"pstate2"])
    expect_equal(as.double(pred.RR2b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)

    # quantile(as.double(pred.RR1a$absRisk.se) - pred.probtrans[,"se2"])
    expect_equal(as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)
    expect_equal(as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)
    expect_equal(as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 5e-3)
    expect_equal(as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 5e-3)
    expect_equal(as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)    
    expect_equal(as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)    
})

## ** With covariates
test_that("predict.CSC (covariates): compare to mstate",{
    for(iX in 0:1){
        newdata <- data.frame(X1 = iX, X2 = 0, X16 = 0)
        newdata.L <- data.frame(X1.1 = c(iX, 0), X1.2 = c(0, iX),
                                X2.1 = c(0, 0), X2.2 = c(0, 0),
                                X16.1 = c(0, 0), X16.2 = c(0, 0),
                                trans = c(1, 2), strata = c(1, 2))
        # mstate
        e.coxph <- coxph(Surv(time, status) ~ X1.1 + X1.2 + X2.1 + X2.2 + X16.1 + X16.2 + strata(trans),
                         data = dL.exp)
        pred.msfit <- msfit(e.coxph, newdata = newdata.L, trans = tmat)
        suppressWarnings(
            pred.probtrans <- probtrans(pred.msfit,0)[[1]]
        )
        ## riskRegression
        CSC.RR1 <- CSC(Hist(time,event)~X1+X2+X16, data = d, method = "breslow")
        pred.RR1a <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
        pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
        pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
        pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
        CSC.RR2 <- CSC(Hist(time,event)~X1+X2+X16, data = d, surv.type = "survival", method = "breslow")
        pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
        pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
        expect_equal(as.double(pred.RR1a$absRisk),pred.probtrans[,"pstate2"])
        expect_equal(as.double(pred.RR1b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)
        expect_equal(as.double(pred.RR1c$absRisk),pred.probtrans[,"pstate3"])
        expect_equal(as.double(pred.RR1d$absRisk),pred.probtrans[,"pstate3"], tol = 1e-2)
        expect_equal(as.double(pred.RR2a$absRisk),pred.probtrans[,"pstate2"], tol = 1e-2)
        expect_equal(as.double(pred.RR2b$absRisk),pred.probtrans[,"pstate2"], tol = 1e-1)
        #if(iX==0){
        # quantile(as.double(pred.RR1a$absRisk.se) - pred.probtrans[,"se2"])
        # quantile(as.double(pred.RR1c$absRisk.se) - pred.probtrans[,"se3"])
        expect_equal(as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        expect_equal(as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        if(iX==0){
            expect_equal(as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 1e-2)
            expect_equal(as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 1e-2)
        }
        expect_equal(as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)    
        expect_equal(as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        ## }else{
        ##     expect_equal(as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)
        ##     expect_equal(as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)
        ##     expect_equal(as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 1e-1)
        ##     expect_equal(as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 1e-1)
        ##     expect_equal(as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)    
        ##     expect_equal(as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)
        ## }
    }
})

## ** Strata 
test_that("predict.CSC (strata): compare to mstate",{

    newdata <- data.frame(X1 = 0, X2 = 0, X16 = 0)
    newdata.L <- data.frame(X1.1 = c(0, 0), X1.2 = c(0, 0),
                            X2.1 = c(0, 0), X2.2 = c(0, 0),
                            X16.1 = c(0, 0), X16.2 = c(0, 0),
                            trans = c(1, 2), strata = c(1, 2))

    e.coxph <- coxph(Surv(time, status) ~ X1.1 + X1.2 + X2.1 + X2.2 + X16.1 + X16.2 + strata(trans),
                     data = dL.exp)
    
    pred.msfit <- msfit(e.coxph, newdata = newdata.L, trans = tmat)
    suppressWarnings(
        pred.probtrans <- probtrans(pred.msfit,0)[[1]]
    )

    ## riskRegression no strata
    CSC.RR1 <- CSC(Hist(time,event)~X1+X2+X16, data = d, method = "breslow")
    pred.RR1a <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
    
    pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
    
    CSC.RR2 <- CSC(Hist(time,event)~X1+X2+X16, data = d, surv.type = "survival", method = "breslow")
    pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
    pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)

    ## riskRegression strata
    d2 <- rbind(cbind(d,grp=1),cbind(d,grp=2))
    CSC.RR1_strata <- CSC(Hist(time,event)~X1+X2+X16+strata(grp), data = d2, method = "breslow")
    pred.RR1a_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1b_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)
    
    pred.RR1c_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)

    pred.RR1d_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE)

    expect_equal(pred.RR1a_strata$absRisk,pred.RR1a$absRisk)
    expect_equal(pred.RR1b_strata$absRisk,pred.RR1b$absRisk)
    expect_equal(pred.RR1c_strata$absRisk,pred.RR1c$absRisk)
    expect_equal(pred.RR1d_strata$absRisk,pred.RR1d$absRisk)

    CSC.RR2_strata<- CSC(Hist(time,event)~X1+X2+X16+strata(grp), data = d2, surv.type = "survival", method = "breslow")
    pred.RR2a_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
    pred.RR2b_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                                keep.newdata = FALSE, se = TRUE, product.limit = FALSE)

    expect_equal(pred.RR2a_strata$absRisk,pred.RR2a$absRisk)
    expect_equal(pred.RR2b_strata$absRisk,pred.RR2b$absRisk)


})
## ** Vs. fixed values
set.seed(10)
n <- 300
df.S <- prodlim::SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("[predictCSC] vs. fixed values",{
    CSC.S <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = "efron", fitter = "coxph")

    ## cause 1
    Event1.S <- predict(CSC.S, newdata = df.S[1:5,], times = seqTime, cause = 1,
                        se = TRUE, keep.newdata = FALSE)

    ## butils::object2script(as.data.table(Event1.S), digit = 6)
    GS <- data.table("observation" = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
           "times" = c( 0.14,  0.14,  0.14,  0.14,  0.14,  0.36,  0.36,  0.36,  0.36,  0.36,  0.56,  0.56,  0.56,  0.56,  0.56,  0.94,  0.94,  0.94,  0.94,  0.94,  3.23,  3.23,  3.23,  3.23,  3.23,  4.12,  4.12,  4.12,  4.12,  4.12,  6.75,  6.75,  6.75,  6.75,  6.75,  9.19,  9.19,  9.19,  9.19,  9.19, 15.90, 15.90, 15.90, 15.90, 15.90, 16.90, 16.90, 16.90, 16.90, 16.90), 
           "strata" = c("X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2"), 
           "absRisk" = c(0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.010614, 0.008321, 0.000000, 0.005059, 0.033321, 0.010614, 0.008321, 0.000000, 0.005059, 0.033321, 0.187045, 0.149948, 0.024562, 0.094195, 0.208568, 0.285231, 0.232052, 0.047504, 0.148967, 0.320621, 0.532564, 0.455475, 0.127580, 0.315670, 0.471993, 0.592726, 0.516242, 0.341738, 0.368635, 0.694506, NA, NA, 0.804087, NA, NA, NA, NA, NA, NA, NA), 
           "absRisk.se" = c(0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.010671, 0.008392, 0.000000, 0.005145, 0.021735, 0.010671, 0.008392, 0.000000, 0.005145, 0.021735, 0.045082, 0.037994, 0.013808, 0.026548, 0.060483, 0.053248, 0.046578, 0.021627, 0.034699, 0.079498, 0.075317, 0.072541, 0.047915, 0.062280, 0.085417, 0.076575, 0.076338, 0.087346, 0.069072, 0.159006, NA, NA, 0.145491, NA, NA, NA, NA, NA, NA, NA))
                         
    expect_equal(as.data.table(Event1.S)$absRisk, GS$absRisk, tolerance = 1e-5)
    expect_equal(as.data.table(Event1.S)$absRisk.se, GS$absRisk.se, tolerance = 1e-5)
    
    ## cause 2
    Event2.S <- predict(CSC.S, newdata = df.S[1:5,], times = seqTime, cause = 2,
                        se = TRUE, keep.newdata = FALSE)
    ## butils::object2script(as.data.table(Event2.S), digit = 6)
    GS <- data.table("observation" = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
                     "times" = c( 0.14,  0.14,  0.14,  0.14,  0.14,  0.36,  0.36,  0.36,  0.36,  0.36,  0.56,  0.56,  0.56,  0.56,  0.56,  0.94,  0.94,  0.94,  0.94,  0.94,  3.23,  3.23,  3.23,  3.23,  3.23,  4.12,  4.12,  4.12,  4.12,  4.12,  6.75,  6.75,  6.75,  6.75,  6.75,  9.19,  9.19,  9.19,  9.19,  9.19, 15.90, 15.90, 15.90, 15.90, 15.90, 16.90, 16.90, 16.90, 16.90, 16.90), 
                     "strata" = c("X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2", "X1=1 X3=1", "X1=1 X3=1", "X1=1 X3=2", "X1=1 X3=1", "X1=0 X3=2"), 
                     "absRisk" = c(0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015510, 0.015342, 0.028621, 0.014982, 0.000000, 0.114820, 0.115596, 0.028621, 0.115825, 0.000000, 0.133554, 0.135236, 0.065157, 0.136670, 0.000000, 0.258482, 0.280051, 0.201121, 0.314747, 0.185563, 0.301218, 0.336703, 0.368645, 0.398416, 0.185563, NA, NA, 0.368645, NA, NA, NA, NA, NA, NA, NA), 
                     "absRisk.se" = c(0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015393, 0.015250, 0.027535, 0.014981, 0.000000, 0.040650, 0.041823, 0.027535, 0.044450, 0.000000, 0.043741, 0.045308, 0.044233, 0.048753, 0.000000, 0.058568, 0.064056, 0.082396, 0.076642, 0.085257, 0.064560, 0.071313, 0.108398, 0.087042, 0.085257, NA, NA, 0.108398, NA, NA, NA, NA, NA, NA, NA))

    expect_equal(as.data.table(Event2.S)$absRisk, GS$absRisk, tolerance = 1e-5)
    expect_equal(as.data.table(Event2.S)$absRisk.se, GS$absRisk.se, tolerance = 1e-5)

})

## ** Melanoma
data(Melanoma, package = "riskRegression")
Melanoma$index <- 1:NROW(Melanoma)

## CSC
cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+sex,
                          Hist(time,status)~age+sex),
             data=Melanoma)
cfit2 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+sex,
                          Hist(time,status)~age+sex),
             data=Melanoma, surv.type = "survival")

## mstrate
Melanoma$event1 <- as.numeric(Melanoma$event == "death.malignant.melanoma")
Melanoma$event2 <- as.numeric(Melanoma$event == "death.other.causes")

MelanomaL <- msprep(time = c(NA, "time", "time"),
                    status = c(NA,"event1", "event2"),
                    data = Melanoma, keep = c("age","logthick","epicel","sex","index"),
                    trans = tmat)
MelanomaL.exp <- expand.covs(MelanomaL,  c("age","logthick","epicel","sex"))

Melanoma.coxph <- coxph(Surv(time, status) ~ age.1 + age.2 + logthick.1 + epicelpresent.1 + sexMale.1 + sexMale.2 + strata(trans),
                        data = MelanomaL.exp)

newdata <- Melanoma[1,,drop=FALSE]
newdata.exp <- MelanomaL.exp[MelanomaL.exp$index %in% newdata$index,,drop=FALSE]
newdata.exp$strata <- newdata.exp$trans

pred.msfit <- msfit(Melanoma.coxph, newdata = newdata.exp, trans = tmat, variance = FALSE)
suppressWarnings(
    pred.probtrans <- probtrans(pred.msfit,0,variance = FALSE)[[1]]
)

## check
test_that("predict.CSC (Melanoma): compare to mstate",{

    pred.RR1 <- predict(cfit1, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = TRUE) 
    expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR1$absRisk))

    pred.RR2 <- predict(cfit1, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
    expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR2$absRisk), tol = 5e-3)

    pred.RR3 <- predict(cfit2, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
    expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR3$absRisk), tol = 1e-1)
})

## * [predictCSC] Conditional absolute risk
cat("[predictCSC] Conditional absolute risk \n")

## ** data and model
set.seed(10)
d <- prodlim::SimCompRisk(1e2)
d$time <- round(d$time,1)
ttt <- sample(x = unique(sort(d$time)), size = 10)
d2 <- prodlim::SimCompRisk(1e2)

CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow")

## ** before the first event
test_that("[predictCSC]: Conditional CIF identical to CIF before first event", {
  pred <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, landmark = min(d$time)-1e-5)
  expect_equal(pred, predC)

  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt)
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, landmark = min(d$time)-1e-5)
  expect_equal(pred, predC)
})

## ** after the last observation
test_that("[predictCSC]: Conditional CIF is NA after the last observation", {
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, landmark = max(d$time)+1)
  expect_equal(all(is.na(predC$absRisk)), TRUE)
  
  predC <- predict(CSC.fit, newdata = d2, cause = 2, times = ttt, landmark = max(d$time)+1)
  expect_equal(all(is.na(predC$absRisk)), TRUE)
  
  t0 <- mean(range(d$time))
  ttt0 <- c(t0,ttt)
  predC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt0, landmark = t0)
  expect_equal(all(is.na(predC$absRisk[,ttt0<t0])), TRUE)
  expect_equal(all(!is.na(predC$absRisk[,ttt0>=t0])), TRUE)
})

## ** at the last event
test_that("[predictCSC] Value of the conditional CIF | at the last event", {
    tau <-  max(d[d$cause==2,"time"])
    
    predC_auto <- predict(CSC.fit, newdata = d2[1:5,], cause = 2, times = tau, landmark = tau, product.limit = FALSE)
    pred <- as.data.table(predictCox(CSC.fit$models[[2]], times = tau, newdata = d2[1:5,], type = "hazard"))
    expect_equal(as.double(predC_auto$absRisk),pred$hazard)

    predC_auto <- predict(CSC.fit, newdata = d2[1:5,], cause = 2, times = tau, landmark = tau, product.limit = TRUE)
    pred <- as.data.table(predictCox(CSC.fit$models[[2]], times = tau, newdata = d2[1:5,], type = "hazard"))
    expect_equal(as.double(predC_auto$absRisk),pred$hazard)
})


## ** vs. manual calculation
test_that("[predictCSC] conditional CIF vs. manual calculation", {
  sttt <- sort(c(0,ttt))
  indexT0 <- 5
    
  # product.limit = FALSE
  cumH1 <- predictCox(CSC.fit$models$`Cause 1`, newdata = d2, times = sttt[indexT0]-1e-6)[["cumhazard"]]
  cumH2 <- predictCox(CSC.fit$models$`Cause 2`, newdata = d2, times = sttt[indexT0]-1e-6)[["cumhazard"]]
  Sall <- exp(-cumH1-cumH2)
  
  predRef <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt[indexT0]-1e-6, product.limit = FALSE)
  
  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, product.limit = FALSE)
  predC_manuel <- (pred$absRisk-as.double(predRef$absRisk))/as.double(Sall)
  predC_manuel[,seq(1,indexT0-1)] <- NA
  
  predC_auto <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, landmark = sttt[indexT0], product.limit = FALSE)
  expect_equal(predC_auto$absRisk,predC_manuel)
  # predC_auto$absRisk - predC_manuel

  # product.limit = TRUE
  h1 <- predictCox(CSC.fit$models$`Cause 1`, newdata = d2, times = CSC.fit$eventTimes, type = "hazard")[["hazard"]]
  h2 <- predictCox(CSC.fit$models$`Cause 2`, newdata = d2, times = CSC.fit$eventTimes, type = "hazard")[["hazard"]]
  Sall <- apply(1-h1-h2,1, function(x){
      c(1,cumprod(x))[sindex(jump.times = CSC.fit$eventTimes, eval.times = sttt[indexT0]-1e-6)+1]
  })
  predRef <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt[indexT0]-1e-6, product.limit = TRUE)
  
  pred <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, product.limit = TRUE)
  predC_manuel <- sweep(pred$absRisk-as.double(predRef$absRisk), MARGIN = 1, FUN ="/", STATS = as.double(Sall))
  predC_manuel[,seq(1,indexT0-1)] <- NA
  
  predC_auto <- predict(CSC.fit, newdata = d2, cause = 2, times = sttt, landmark = sttt[indexT0], product.limit = TRUE)
  expect_equal(predC_auto$absRisk,predC_manuel)
  # predC_auto$absRisk-predC_manuel
})


## ** vs. fixed values
set.seed(10)
d <- sampleData(3e2, outcome = "competing.risks")
d$time <- round(d$time,2)
ttt <- sample(x = unique(sort(d$time)), size = 10) 
d2 <- d
times <- sort(c(0,d$time))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("[predictCSC] conditional CIF vs. fixed values",{
  CSC.fitS <- CSC(Hist(time,event)~ strata(X1) + X5 + strata(X3) + X7 +X2,data=d, method = "breslow", surv.type = "survival")

  p1 <- predict(CSC.fitS, newdata = d2[1:5,], times = seqTime, landmark = 1, cause = 1,
                se = FALSE, keep.newdata = FALSE)

    ## butils::object2script(as.data.table(p1), digit = 6)

  GS <- data.table("observation" = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
                   "times" = c( 0.14,  0.14,  0.14,  0.14,  0.14,  0.36,  0.36,  0.36,  0.36,  0.36,  0.56,  0.56,  0.56,  0.56,  0.56,  0.94,  0.94,  0.94,  0.94,  0.94,  3.23,  3.23,  3.23,  3.23,  3.23,  4.12,  4.12,  4.12,  4.12,  4.12,  6.75,  6.75,  6.75,  6.75,  6.75,  9.19,  9.19,  9.19,  9.19,  9.19, 15.90, 15.90, 15.90, 15.90, 15.90, 16.90, 16.90, 16.90, 16.90, 16.90), 
                   "strata" = c("X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0", "X1=0 X3=0"), 
                   "absRisk" = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.193716, 0.321904, 0.323271, 0.158113, 0.148124, 0.258450, 0.422776, 0.423841, 0.212179, 0.197920, 0.344406, 0.547189, 0.546861, 0.285912, 0.264501, 0.405789, 0.628009, 0.625925, 0.340374, 0.312468, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

    expect_equal(as.data.table(p1)$absRisk, GS$absRisk, tolerance = 1e-5)
    expect_equal(as.data.table(p1)$absRisk.se, GS$absRisk.se, tolerance = 1e-5)

  
  expect_error(predict(CSC.fitS, newdata = d2[1:10,], times = seqTime, landmark = 1,
                       cause = 1, se = TRUE))
  ## Error in predict.CauseSpecificCox(CSC.fitS, newdata = d2[1:10, ], times = seqTime,  : 
  ## standard error for the conditional survival not implemented 
})


## * [predictCSC] Argument diag
cat("[predictCox] Argument \'diag\' \n")
set.seed(10)
dt <- sampleData(75, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]

test_that("[predictCox] diag no strata", {
    e.CSC <- CSC(Hist(time, event) ~ X1*X6, data = dt)

    GS <- predict(e.CSC, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, average.iid = TRUE, cause = 1)
    test <- predict(e.CSC, newdata = dt, times = dt$time,
                    se = FALSE, iid = TRUE, average.iid = TRUE, diag = TRUE, cause = 1)
    test2 <- predict(e.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE, cause = 1)
    
    ## estimates
    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$absRisk), as.double(test$absRisk))

    ## iid
    GS.iid.diag <- do.call(rbind,lapply(1:NROW(dt),
                                        function(iN){GS$absRisk.iid[iN,iN,]}))
    expect_equal(GS.iid.diag, test$absRisk.iid[,1,])

    ## average.iid
    expect_equal(colMeans(GS.iid.diag), test2$absRisk.average.iid[,1])    
    expect_equal(test$absRisk.average.iid, test2$absRisk.average.iid)
    range(test$absRisk.average.iid - colMeans(GS.iid.diag))

    test2 <- predict(e.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE, cause = 1)
    dim(test2$absRisk.average.iid)

})

test_that("[predictCox] diag strata", {
    eS.CSC <- CSC(Hist(time, event) ~ strata(X1) + X6, data = dt)

    GS <- predict(eS.CSC, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, cause = 1)
    test <- predict(eS.CSC, newdata = dt, times = dt$time, se = FALSE, iid = TRUE, diag = TRUE, cause = 1)

    expect_equal(dt$time, as.double(test$time))
    expect_equal(diag(GS$absRisk), as.double(test$absRisk))

    GS.iid.diag <- do.call(rbind,lapply(1:NROW(dt),
                                        function(iN){GS$absRisk.iid[iN,iN,]}))
    expect_equal(GS.iid.diag, test$absRisk.iid[,1,])

})

## * [predictCSC] Average iid
cat("[predictCSC] Average iid \n")
set.seed(10)
d <- sampleData(100, outcome = "competing.risks")

obs.firstEvent <- which.min(d$time)
strata.firstEvent <- 3

seqTime <- unique(c(0,min(d$time),d$time[sample.int(n = 100, size = 10)],1e8))
index.firstEvent <- 2

## ** non parametric

test_that("iid average - non parametric (hazard)", {
    for(iType in c("hazard","survival")){ ## iType <- "hazard"
        m.CSC <- CSC(Hist(time, event) ~ strata(X1) + strata(X2),
                     data = d, surv.type = iType)

        res1 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, iid = TRUE, average.iid = TRUE)

        res2 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, average.iid = TRUE)

        expect_equal(res1$absRisk.average.iid,res2$absRisk.average.iid)
        expect_true(all(res1$absRisk.average.iid[,1]==0))
        expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))

        expect_equal(t(apply(res1$absRisk.iid,2:3,mean)),
                     res2$absRisk.average.iid)
        
        ## compare to fixed value    
        ## d[time==min(time),]
        ## levels(predictCox(m.CSC$models[[1]])$strata)
        expect_equal(res1$absRisk.iid[obs.firstEvent, index.firstEvent,],
                     iidCox(m.CSC$models[[1]], return.object = FALSE)$IFhazard[[strata.firstEvent]][,1])
    }
})

## ** semi parametric
test_that("iid average - semi parametric", {
    for(iType in c("hazard","survival")){ ## iType <- "hazard"
        m.CSC <- CSC(Hist(time, event) ~ X1*X6 + strata(X2), data = d, surv.type = iType)
        
        res1 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, iid = TRUE, average.iid = TRUE)
        res2 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, average.iid = TRUE)
        GS3 <- res1$absRisk.average.iid
        GS23 <- res1$absRisk.average.iid
        
        expect_equal(res1$absRisk.average.iid,res2$absRisk.average.iid)
        expect_true(all(res1$absRisk.average.iid[,1]==0))
        expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))

        expect_equal(t(apply(res1$absRisk.iid,2:3,mean)),
                     res2$absRisk.average.iid)

    }

})

## * [predictCSC] Miscelaneous
## ** Confidence bands vs. timereg
cat("[predictCSC] Confidence band vs. timereg \n")

## ** Data
set.seed(10)
dt <- sampleData(1e2, outcome = "competing.risks")
newdata <- dt[1:10,]

## ** Model
e.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = dt)
vec.times <- e.CSC$eventTimes

e.timereg <- comp.risk(Event(time, event) ~ const(X1) + const(X2),
                       model = "rcif",
                       data = dt, cause = 1)
coef(e.timereg)
coef(e.CSC)

## ** Compute confidence bands
if(FALSE){
    resTimereg <- predict.timereg(e.timereg,
                                  newdata = newdata[1,,drop=FALSE],
                                  times = vec.times,
                                  resample.iid = 1,
                                  n.sim = nsim.band)
    resTimereg$P1
}
predRR <- predict(e.CSC,
                  newdata = newdata,
                  times = vec.times-1e-5,
                  se = TRUE,
                  band = TRUE,
                  nsim.band = nsim.band,
                  cause = 1)
## debuging

## cumHazard.coxph <- predictCox(fit.coxph)$cumhazard
## iid.coxph <- iidCox(fit.coxph)
## iid.lambda <- iid.coxph$ICcumhazard[[1]][1,]

## X.design <- model.matrix(formula(fit.coxph),newdata[i,,drop=FALSE])[,-1]
## eLP <- exp(X.design %*% cbind(coef(fit.coxph)))
## Xiid.beta <- X.design %*% iid.coxph$IFbeta[1,]

## term1 <- as.numeric(eLP * iid.lambda)
## term2 <- as.numeric(eLP * cumHazard.coxph * Xiid.beta)

## ls.args <- list(delta = res$cumhazard.iid[i,,],
##                 nObs = 1,
##                 nt = length(times),
##                 n = NROW(d),
##                 mpt = n.sim*1,
##                 nSims = n.sim)
## ls.args$se <-  apply(ls.args$delta^2, 1, sum)^0.5

## mpt <- .C("confBandBasePredict", delta = as.double(delta), 
##           nObs = as.integer(nobs), nt = as.integer(nt), n = as.integer(n), 
##           se = as.double(se), mpt = double(n.sim * nobs), nSims = as.integer(n.sim), 
##           PACKAGE = "timereg")$mpt

## ** Order of the prediction times
cat("[predictCSC] order of the prediction times \n")
data(Melanoma, package = "riskRegression")
times2 <- sort(c(0,0.9*min(Melanoma$time),Melanoma$time[5],max(Melanoma$time)*1.1))
newOrder <- sample.int(length(times2),length(times2),replace = FALSE)

test_that("Prediction with CSC - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick, data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1, keep.times = FALSE)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1, keep.times = FALSE)
  expect_equal(predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
})

test_that("Prediction with CSC (strata) - sorted vs. unsorted times",{
  fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion), data = Melanoma)
  predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1)
  predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1)
  expect_equal(predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
})


## ** Dependence on other arguments
cat("[predictCSC] dependence on other arguments \n")
test_that("[predictCSC] output of average.iid should not depend on other arguments", {
    set.seed(10)
    d <- sampleData(70,outcome="competing.risks")
    d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]

    fit <- CSC(Hist(time,event)~X1 + strata(X2) + X6,
               data=d)

    out1 <- predict(fit, newdata = d[1:5], times = 1:3, average.iid = TRUE, cause = 1)
    out2 <- predict(fit, newdata = d[1:5], times = 1:3, se = TRUE, average.iid = TRUE, cause = 1)

    test_that("output of average.iid should not depend on other arguments", {
        expect_equal(out1$survival.average.iid,out2$survival.average.iid, tol = 1e-8)
    })
})

## ** Prediction when iid is stored
cat("[predictCSC] prediction after storing iid \n")
data(Melanoma, package = "riskRegression")
cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+strata(sex),
                          Hist(time,status)~age+strata(sex)),
             data=Melanoma)
test_that("[predictCSC] prediction after storing iid", {
    GS <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
                   times=4,se=TRUE,band=TRUE)
    cfit1 <- iidCox(cfit1)
    res <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
                   times=4,se=TRUE,band=TRUE)
    
    expect_equal(GS,res)
})

## ** Prediction first/NA/negative/last event
cat("[predictCSC] prediction first/NA/negative/last event \n")

## *** no strata
data(Melanoma, package = "riskRegression")
times1 <- unique(Melanoma$time)
times2 <- c(0,0.9*min(times1),times1*1.1)
dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]

fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma, fitter = "cph")

test_that("[predictCSC] (no strata): NA after last event",{
  test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.vector(is.na(prediction$absRisk)), c(FALSE, FALSE, TRUE))
})

test_that("[predictCSC] (no strata): no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction$absRisk), 0)
})

test_that("[predictCSC] (no strata): beyond last event time",{
    test.times.1 <- c(10,3000,5000)
    newd <- Melanoma
    newd <- Melanoma[1,,drop=FALSE]
    prediction.1 <- predict(fit.CSC, times = test.times.1, newdata = newd, cause = 1)
    test.times.2 <- c(10,300,5000)
    prediction.2 <- predict(fit.CSC, times = test.times.2, newdata = newd, cause = 1)
    expect_equal(prediction.1$absRisk[,3],prediction.2$absRisk[,3])
})

test_that("[predictCSC]: Deal with negative time points",{
    expect_equal(unname(predict(fit.CSC, times = -1, newdata = newd, cause = 1)$absRisk),
                 matrix(0,nrow = nrow(newd), ncol = 1))
})

test_that("[predictCSC]:Deal with NA in times",{
  expect_error(predict(fit.CSC, times = c(test.times,NA), newdata = newd, cause = 1))
})

## *** strata
fit.coxph <- coxph(Surv(time,status==1) ~ thick + strata(invasion) + strata(ici), data = Melanoma,
                   x = TRUE, y = TRUE)
fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion) + strat(ici), data = Melanoma, fitter = "cph")

data.test <- data.table(Melanoma)[, .SD[1], by = c("invasion", "ici")]
setkeyv(data.test, c("invasion","ici"))

# identify the last event time for each strata
epsilon <- min(diff(unique(fit.coxph$y[,"time"])))/10
pred.coxph <- predictCox(fit.coxph, keep.strata = TRUE, keep.times = TRUE)
baseHazStrata <- as.data.table(pred.coxph[c("times","hazard","cumhazard","strata","survival")])
dt.times <- baseHazStrata[, .(beforeLastTime = times[.N]-epsilon,
                              LastTime = times[.N],
                              afterLastTime = times[.N]+epsilon),
                          by = strata]

test_that("[predictCSC](strata): NA after last event",{
  for(Ttempo in 1:nrow(dt.times)){
    test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
    
    prediction <- predict(fit.CSC, times = test.times, newdata = data.test, cause = 1)
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
  }
})

test_that("[predictCSC](strata)  - no event before prediction time",{
  test.times <- min(Melanoma$time)-1e-5
  
  prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
  expect_equal(as.double(prediction$absRisk), 0)
})

## ** Argument store.iid
cat("[predictCSC] Argument \'store.iid\' \n")
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")
setkey(d,time)

m.CSC <- CSC(Hist(time, event) ~ X1+X6, data = d)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d
head(newdata[,.(time,event,X1,X6)])

test_that("[predictCSC]: iid minimal - no strata", {
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata[1],
                    cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata[1],
                    cause = 1,
                    store.iid = "minimal", average.iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata[1],
                    cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)

    
    expect_equal(res1$absRisk.se,res3$absRisk.se)
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    expect_equal(res2$absRisk.average.iid, t(apply(res3$absRisk.iid,2:3,mean)))
})


test_that("[predictCSC]: iid minimal - strata", {
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "minimal", average.iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$absRisk.se,res3$absRisk.se)
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    expect_equal(res2$absRisk.average.iid, t(apply(res3$absRisk.iid,2:3,mean)))
})
## * [predictCSC] Possible issue (estimated absolute risk over 1)
## this section does not perform any tests
## but show an example where the estimated absolute risk
## is over 1, probably because of the small sample size.
## I don't know if this is an issue.
if(FALSE){
    set.seed(5)
    d <- sampleData(80,outcome="comp")
    d[event==0,event:=1] # remove censoring
    ttt <- sort(unique(d[X1==1,time]))

    CSC.fit.s <- CSC(list(Hist(time,event)~ strata(X1)+X2+X9,
                          Hist(time,event)~ X2+strata(X4)+X8+X7),data=d)
    predict(CSC.fit.s,newdata=d[X1==1][1],cause=1,times=ttt,se=1L) # NA values due to absolute risk over 1
    predict(CSC.fit.s,newdata=d[X1==1][1],product.limit = FALSE,cause=1,times=ttt,se=1L) # NA values due to absolute risk over 1

    #### investigate how come we get absRisk > 1
    # since absRisk = int hazard1 * survival
    # it is possible if survival > 1
    #                or hazard(1) > 1
    
    ## restrict to the first strata and simplify the model
    ## to see if we can get an hazard > 1
    dt <- d[X1==1]
    setkeyv(dt, "time")
    m.coxph <- coxph(Surv(time,event>=1)~X2, data = dt, x = TRUE, ties = "breslow")

    # compute baseline hazard
    as.data.table(predictCox(m.coxph, centered = TRUE, type = "hazard")) # automatic
    rev(1/cumsum(rev(eXB))) # manual

    # normally they are normalized such that they are at most one:
    rev(rev(eXB)/cumsum(rev(eXB)))

    # but this does not work for new observations
    rev(unique(eXB)[1]/cumsum(rev(eXB))) # ok
    rev(unique(eXB)[2]/cumsum(rev(eXB))) # no hazard over 1

    # this also creates a problem when computing the suvival using the product limit estimator 
}


#----------------------------------------------------------------------
### test-predictCSC.R ends here
