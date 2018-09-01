### test-predictCSC_vs_mstate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: sep  1 2018 (11:24) 
##           By: Brice Ozenne
##     Update #: 193
#----------------------------------------------------------------------
## 
### Commentary: 
## Compare predictCSC with mstate
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

## * 1- Compare predict.CSC and mstate on small examples
cat("[predictCSS] - small example \n")
## ** data
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

## ** 1a- no risk factor no strata
test_that("predict.CSC/predictCox (1a,df1) - compare to manual estimation",{
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

test_that("predict.CSC/predictCox (1a,df2) - compare to manual estimation",{
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

test_that("predict.CSC (1a) - compare to mstate",{
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


test_that("predict.CSC(1a) - add up to one",{
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
# }}}

## ** 1b- no risk factor strata
test_that("predict.CSC(1b) - compare to manual estimation",{

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

## ** 1c- risk factor no strata
test_that("predict.CSC(1c,df1) - compare to manual estimation",{
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

test_that("predict.CSC(1c,df2) - compare to manual estimation",{
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


test_that("predict.CSC(1c) - compare to mstate",{
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

test_that("predict.CSC(1c) - add up to one",{
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

## * 2- compare predict.CSC and mstate on simulated data
cat("[predictCSS] - simulated data \n")
## ** data
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]

d <- d[event!=0]

d[, event1 := as.numeric(event == 1)]
d[, event2 := as.numeric(event == 2)]

## ** prepare mstrate
tmat <- trans.comprisk(2, names = c("0", "1", "2"))

dL <- msprep(time = c(NA, "time", "time"),
             status = c(NA,"event1", "event2"),
             data = d, keep = c("X1","X2","X16","Xcat2"),
             trans = tmat)
dL.exp <- expand.covs(dL,  c("X1","X2","X16","Xcat2"))

## ** 2a - no covariates
test_that("predict.CSC(2a) - compare to mstate",{
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

## ** 2b - with covariates
test_that("predict.CSC(2b) - compare to mstate",{
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

## ** 2c - strata 
test_that("predict.CSC(2a) - compare to mstate",{

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

## * 3- compare predict.CSC and mstate on Melanoma
cat("[predictCSS] - Melanoma \n")
## ** Data
data(Melanoma, package = "riskRegression")
# sum(is.na(Melanoma))
Melanoma$index <- 1:NROW(Melanoma)

cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+sex,
                          Hist(time,status)~age+sex),
             data=Melanoma)
cfit2 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+sex,
                          Hist(time,status)~age+sex),
             data=Melanoma, surv.type = "survival")


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

pred.RR1 <- predict(cfit1, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = TRUE) 
expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR1$absRisk))

pred.RR2 <- predict(cfit1, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR2$absRisk), tol = 5e-3)

pred.RR3 <- predict(cfit2, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
expect_equal(as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR3$absRisk), tol = 1e-1)

## * 4- test predict.CSC with strata
cat("[predictCSS] - strata \n")
## ** CIF
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("predictCSC with strata",{
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

## ** conditional CIF
set.seed(10)
d <- sampleData(3e2, outcome = "competing.risks")
d$time <- round(d$time,2)
ttt <- sample(x = unique(sort(d$time)), size = 10) 
d2 <- d
times <- sort(c(0,d$time))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("conditional predictCSC with strata",{
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


## * 5- Fast SE
cat("[predictCSC] - fast SE \n")
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")
setkey(d,time)

m.CSC <- CSC(Hist(time, event) ~ X1+X6, data = d)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

test_that("iid minimal - no strata", {
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

test_that("iid minimal - strata", {
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

## * 6- Average iid
cat("[predictCSC] - average iid \n")
set.seed(10)
d <- sampleData(100, outcome = "competing.risks")

obs.firstEvent <- which.min(d$time)
strata.firstEvent <- 3

seqTime <- unique(c(0,min(d$time),d$time[sample.int(n = 100, size = 10)],1e8))
index.firstEvent <- 2

## ** non parametric
m.CSC <- CSC(Hist(time, event) ~ strata(X1) + strata(X2),
             data = d, surv.type = "hazard")

res2 <- predict(m.CSC, times = seqTime, newdata = d,
                cause = 1, average.iid = TRUE)

test_that("iid average - non parametric (hazard)", {
    
    res1 <- predict(m.CSC, times = seqTime, newdata = d,
                    cause = 1, iid = TRUE, average.iid = TRUE)

    expect_equal(res1$absRisk.average.iid,res2$absRisk.average.iid)
    expect_true(all(res1$absRisk.average.iid[,1]==0))
    expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))

    ## compare to fixed value    
    ## d[time==min(time),]
    ## levels(predictCox(m.CSC$models[[1]])$strata)
    expect_equal(res1$absRisk.iid[obs.firstEvent, index.firstEvent,],
                 iidCox(m.CSC$models[[1]])$IFhazard[[strata.firstEvent]][,1])

})

m.CSC <- CSC(Hist(time, event) ~ strata(X1) + strata(X2),
             data = d, surv.type = "survival")

res2 <- predict(m.CSC, times = seqTime, newdata = d,
                cause = 1, average.iid = TRUE)

test_that("iid average - non parametric (survival)", {

    res1 <- predict(m.CSC, times = seqTime, newdata = d,
                    cause = 1, iid = TRUE, average.iid = TRUE)


    expect_equal(res1$absRisk.average.iid,res2$absRisk.average.iid)
    expect_true(all(res1$absRisk.average.iid[,1]==0))
    expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))

    ## compare to fixed value    
    ## d[time==min(time),]
    ## levels(predictCox(m.CSC$models[[1]])$strata)
    expect_equal(res1$absRisk.iid[obs.firstEvent, index.firstEvent,],
                 iidCox(m.CSC$models[[1]])$IFhazard[[strata.firstEvent]][,1])

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
    }

})



## * 7- Prediction when iid is stored
cat("[predictCSS] - iid \n")
data(Melanoma, package = "riskRegression")
cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+strata(sex),
                          Hist(time,status)~age+strata(sex)),
             data=Melanoma)
cfit1$iid <- lapply(cfit1$model,iidCox)


res <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
               times=4,se=TRUE,band=TRUE)

#----------------------------------------------------------------------
### test-predictCSC_vs_mstate.R ends here
