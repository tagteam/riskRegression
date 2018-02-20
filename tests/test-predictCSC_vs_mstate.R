### test-predictCSC_vs_mstate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: Feb 19 2018 (18:22) 
##           By: Thomas Alexander Gerds
##     Update #: 143
#----------------------------------------------------------------------
## 
### Commentary: 
## Compare predictCSC with mstate
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(testthat)
library(mstate)
tmat <- trans.comprisk(2, names = c("0", "1", "2"))
library(survival)

# {{{ 1- compare predict.CSC and mstate on small examples

## simulatenous events
df1 <- data.frame(time = rep(1:10,2),
                 event = c(rep(1,10),rep(2,10)),
                 X1 = 0:1
                 )
df1$event1 <- as.numeric(df1$event == 1)
df1$event2 <- as.numeric(df1$event == 2)

dfS <- rbind(cbind(df1, grp = 1, X2 = 0),
             cbind(rbind(df1,df1),grp = 2, X2 = 0),
             cbind(df1, grp = 3, X2 = 0)
             )

## distinct events
df2 <- data.frame(time = c(1:10,-0.1+1:10),
                  event = c(rep(1,10),rep(2,10)),
                  X1 = 0:1
                  )
df2$event1 <- as.numeric(df2$event == 1)
df2$event2 <- as.numeric(df2$event == 2)

# {{{ 1a- no risk factor no strata

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

# {{{ 1b- no risk factor strata

test_that("predict.CSC(1b) - compare to manual estimation",{
    # don't work:
    # CSC.exp <- CSC(Hist(time,event) ~ strata(grp), data = dfS)
    # X2 factive variable
    CSC.exp <- CSC(Hist(time,event) ~ strata(grp) + X2, data = dfS, method = "breslow")
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(grp = 1:3, X2 = 0), cause = 1, product.limit = FALSE)  
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(grp = 1:3, X2 = 0), cause = 1, product.limit = TRUE)  

    CSC.prodlim <- CSC(Hist(time,event) ~ strata(grp) + X2, data = dfS, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(grp = 1:3, X2 = 0), cause = 1, product.limit = TRUE)  
    expect_equal(pred.exp2,pred.prodlim)
    
    CSC.prodlimE <- CSC(Hist(time,event) ~ strata(grp) + X2, data = dfS, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(grp = 1:3, X2 = 0), cause = 1, product.limit = TRUE)  

    ## baseline
    lambda1 <- lambda2 <- 1/seq(20,2,by = -2)
    lambda2E.1 <- as.data.table(predictCox(CSC.prodlimE$models[[1]], type = "hazard"))
    lambda2E.2 <- as.data.table(predictCox(CSC.prodlimE$models[[2]], type = "hazard"))
        
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

    survival <- c(1,cumprod(1-lambda2E.2[strata=="grp=1"][["hazard"]]))[1:10]
    expect_equal(as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=1"]),
                 cumsum(lambda1*survival))
    survival <- c(1,cumprod(1-lambda2E.2[strata=="grp=2"][["hazard"]]))[1:10]
    expect_equal(as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=2"]),
                 cumsum(lambda2E.1[strata=="grp=2"][["hazard"]]*survival))    
    expect_equal(as.double(pred.prodlimE$absRisk[1,]),
                 as.double(pred.prodlimE$absRisk[3,]))
})
# }}}

# {{{ 1c- risk factor no strata
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

# }}}

# }}}

# {{{ 2- compare predict.CSC and mstate on simulated data


# {{{ data
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]

d <- d[event!=0]

d[, event1 := as.numeric(event == 1)]
d[, event2 := as.numeric(event == 2)]
# }}}

# {{{ prepare mstrate
tmat <- trans.comprisk(2, names = c("0", "1", "2"))

dL <- msprep(time = c(NA, "time", "time"),
             status = c(NA,"event1", "event2"),
             data = d, keep = c("X1","X2","X16","Xcat2"),
             trans = tmat)
dL.exp <- expand.covs(dL,  c("X1","X2","X16","Xcat2"))
# }}}

# {{{ 2a - no covariates
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
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)

    pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
    
    CSC.RR2 <- CSC(Hist(time,event)~1, data = d, surv.type = "survival", method = "breslow")
    pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                        keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                        keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)


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
# }}}

# {{{ 2b - with covariates
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
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)
        pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
        pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)
        pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
        CSC.RR2 <- CSC(Hist(time,event)~X1+X2+X16, data = d, surv.type = "survival", method = "breslow")
        pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)
        pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                             keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
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

# }}}


# {{{ 2c - strata 
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
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1b <- predict(CSC.RR1, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
    
    pred.RR1c <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1d <- predict(CSC.RR1, newdata, cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
    
    CSC.RR2 <- CSC(Hist(time,event)~X1+X2+X16, data = d, surv.type = "survival", method = "breslow")
    pred.RR2a <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)
    pred.RR2b <- predict(CSC.RR2, newdata, cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)

    ## riskRegression strata
    d2 <- rbind(cbind(d,grp=1),cbind(d,grp=2))
    CSC.RR1_strata <- CSC(Hist(time,event)~X1+X2+X16+strata(grp), data = d2, method = "breslow")
    pred.RR1a_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1b_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)
    
    pred.RR1c_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)

    pred.RR1d_strata <- predict(CSC.RR1_strata, cbind(newdata,grp=1), cause = 2, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)

    expect_equal(pred.RR1a_strata$absRisk,pred.RR1a$absRisk)
    expect_equal(pred.RR1b_strata$absRisk,pred.RR1b$absRisk)
    expect_equal(pred.RR1c_strata$absRisk,pred.RR1c$absRisk)
    expect_equal(pred.RR1d_strata$absRisk,pred.RR1d$absRisk)

    CSC.RR2_strata<- CSC(Hist(time,event)~X1+X2+X16+strata(grp), data = d2, surv.type = "survival", method = "breslow")
    pred.RR2a_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE, log.transform = FALSE)
    pred.RR2b_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                                keep.newdata = FALSE, se = TRUE, product.limit = FALSE, log.transform = FALSE)

    expect_equal(pred.RR2a_strata$absRisk,pred.RR2a$absRisk)
    expect_equal(pred.RR2b_strata$absRisk,pred.RR2b$absRisk)


})

# }}}
# }}}

# {{{ 3- compare predict.CSC and mstate on Melanoma
data(Melanoma)
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

# }}}

# {{{ 4- test predict.CSC with strata

exportRes <- function(pred){        
    mat <- apply(pred, 2, paste0, ",")
    mat[,NCOL(mat)] <- gsub(",","),",mat[,NCOL(mat)])
    mat[,1] <- paste0("c(",mat[,1])
    mat[1,1] <- paste0("rbind(",mat[1,1])
    mat[NROW(mat),NCOL(mat)] <- gsub(",",")",mat[NROW(mat),NCOL(mat)])
    mat <- apply(mat,1,paste, collapse=" ")
    print(as.data.frame(mat), row.names = FALSE, quote = FALSE)
}

# {{{ CIF
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("predictCSC with strata",{
    CSC.S <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = "efron", fitter = "coxph")

    ## cause 1
    Event0.S <- predict(CSC.S, newdata = df.S[1:10,], times = seqTime, cause = 1, se = TRUE, log.transform = FALSE)
    # exportRes(Event0.S$absRisk)
    EventTest.S <- rbind(c(0, 0, 0.0106144145911254, 0.0106144145911254, 0.187045453946268, 0.285231447678766, 0.532564003644496, 0.592725502528842, NA, NA),
                         c(0, 0, 0.00832104772556249, 0.00832104772556249, 0.149947950248539, 0.232051997221862, 0.455475368419988, 0.516242447145683, NA, NA),
                         c(0, 0, 0, 0, 0.0245616111585614, 0.0475040308253248, 0.127579540721806, 0.341737537836591, 0.804087374199052, NA),
                         c(0, 0, 0.00505938518189251, 0.00505938518189251, 0.0941951131187111, 0.148966970722694, 0.315670284337352, 0.368634945834961, NA, NA),
                         c(0, 0, 0.0333211699714726, 0.0333211699714726, 0.208568364123239, 0.320620669066337, 0.471992949590838, 0.694505610551621, NA, NA),
                         c(0, 0, 0.0373150558496205, 0.0373150558496205, 0.230954232455666, 0.352180064158924, 0.511945045355676, 0.735072611361283, NA, NA),
                         c(0, 0, 0, 0, 0.0298096235333887, 0.0574993211690947, 0.152877662070974, 0.39517546530519, 0.85144688504677, NA),
                         c(0, 0, 0, 0, 0.0464226204403023, 0.0810285848692912, 0.200212324558954, 0.280410764048454, NA, NA),
                         c(0, 0, 0.00147567889147631, 0.00147567889147631, 0.0285592182677721, 0.0463338030168861, 0.107882405507319, 0.131153683837835, NA, NA),
                         c(0, 0, 0.00763069515878888, 0.00763069515878888, 0.138447165088275, 0.215222599582239, 0.429105291382415, 0.489238856085115, NA, NA)
                         )
                         
    # exportRes(Event0.S$absRisk.se)
    EventTest.Sse <- rbind(c(0, 0, 0.0106712159100338, 0.0106712159100338, 0.0450817300107038, 0.0532480459954564, 0.0753166210562018, 0.0765745873865227, NA, NA),
                           c(0, 0, 0.00839193402284965, 0.00839193402284965, 0.0379944475244501, 0.0465778404253987, 0.0725406799359907, 0.0763377529996142, NA, NA),
                           c(0, 0, 0, 0, 0.0138081037138213, 0.0216274430722753, 0.0479147789638878, 0.0873460623038645, 0.145491084163842, NA),
                           c(0, 0, 0.00514459679767462, 0.00514459679767462, 0.0265484550995743, 0.0346993581017717, 0.062279517239144, 0.069071517387357, NA, NA),
                           c(0, 0, 0.021735230875603, 0.021735230875603, 0.060482801398461, 0.0794977001947496, 0.0854169134837956, 0.159006165727376, NA, NA),
                           c(0, 0, 0.0242983977238479, 0.0242983977238479, 0.0657251426329389, 0.0845674863329775, 0.0881208514994792, 0.159585617247264, NA, NA),
                           c(0, 0, 0, 0, 0.0165989593694511, 0.0257185232695037, 0.0556061384050834, 0.0952602295190898, 0.160166340233958, NA),
                           c(0, 0, 0, 0, 0.0227165681714911, 0.0303578180154916, 0.0472221634249264, 0.0776650874848463, NA, NA),
                           c(0, 0, 0.00154615285626028, 0.00154615285626028, 0.0105322698341324, 0.0152538305227429, 0.0325610788658487, 0.0385774537048457, NA, NA),
                           c(0, 0, 0.00770540094368359, 0.00770540094368359, 0.0357230071460248, 0.0443332122465385, 0.0711093551906959, 0.0756423921566412, NA, NA))
    
    expect_equal(as.double(Event0.S$absRisk),as.double(EventTest.S), tolerance = 1e-8)
    expect_equal(as.double(Event0.S$absRisk.se),as.double(EventTest.Sse), tolerance = 1e-8)

    ## cause 2
    Event0.S <- predict(CSC.S, newdata = df.S[1:10,], times = seqTime, cause = 2, se = TRUE, log.transform = FALSE)
    # exportRes(Event0.S$absRisk)
    EventTest.S <- rbind(c(0, 0, 0, 0.0155103901004989, 0.114819674310264, 0.133554023880907, 0.258482028753579, 0.30121775499437, NA, NA),
                         c(0, 0, 0, 0.0153421834906167, 0.115596334553523, 0.135235912444512, 0.280051350712365, 0.336702700208684, NA, NA),
                         c(0, 0, 0, 0.0286212020651688, 0.0286212020651688, 0.0651572997890421, 0.201121040112556, 0.368644727533912, 0.368644727533912, NA),
                         c(0, 0, 0, 0.0149823141990889, 0.115825276376062, 0.136670102995648, 0.314747481972403, 0.398416184718663, NA, NA),
                         c(0, 0, 0, 0, 0, 0, 0.185562718904235, 0.185562718904235, NA, NA),
                         c(0, 0, 0, 0, 0, 0, 0.174414830957182, 0.174414830957182, NA, NA),
                         c(0, 0, 0, 0.0289271011335303, 0.0289271011335303, 0.0656370843625444, 0.200372699021718, 0.355456697520017, 0.355456697520017, NA),
                         c(0, 0, 0, 0, 0.0736559703796315, 0.13509415399219, 0.327724288675815, 0.398138259627919, NA, NA),
                         c(0, 0, 0, 0.0140631051282532, 0.112074343769209, 0.133601796979203, 0.351866109834586, 0.476952467684156, NA, NA),
                         c(0, 0, 0, 0.0152808234903303, 0.115754081114025, 0.135661197498255, 0.287027747011357, 0.348624181987888, NA, NA))
    
    # exportRes(Event0.S$absRisk.se)
    EventTest.Sse <-  rbind(c(0, 0, 0, 0.0153926104287061, 0.0406498339336404, 0.0437413700655553, 0.058568277002433, 0.064559589244607, NA, NA),
                            c(0, 0, 0, 0.0152497691155128, 0.041822870900744, 0.0453079920798889, 0.0640558545641551, 0.0713128040165271, NA, NA),
                            c(0, 0, 0, 0.0275346162247003, 0.0275346162247003, 0.044232719129842, 0.0823959408937491, 0.108398160699563, 0.108398160699563, NA),
                            c(0, 0, 0, 0.0149814930865529, 0.0444501861732344, 0.048753408743409, 0.0766416176894157, 0.0870424815241679, NA, NA),
                            c(0, 0, 0, 0, 0, 0, 0.0852567514611162, 0.0852567514611162, NA, NA),
                            c(0, 0, 0, 0, 0, 0, 0.0814413109173167, 0.0814413109173167, NA, NA),
                            c(0, 0, 0, 0.0278920997358414, 0.0278920997358414, 0.0445885438001586, 0.0816345326243422, 0.104724109635943, 0.104724109635943, NA),
                            c(0, 0, 0, 0, 0.0350870280714417, 0.0472595265202325, 0.0747537082009305, 0.0822009076753997, NA, NA),
                            c(0, 0, 0, 0.0144861108401961, 0.0521405476323735, 0.0585393957249475, 0.109847679887477, 0.129934423295299, NA, NA),
                            c(0, 0, 0, 0.015200608407024, 0.042257604261811, 0.0458838711097693, 0.066139740199052, 0.0738888072110129, NA, NA))
    
    expect_equal(as.double(Event0.S$absRisk),as.double(EventTest.S), tolerance = 1e-8)
    expect_equal(as.double(Event0.S$absRisk.se),as.double(EventTest.Sse), tolerance = 1e-8)

})
# }}}

# {{{ conditional CIF

set.seed(10)
d <- sampleData(3e2, outcome = "competing.risks")
d$time <- round(d$time,2)
ttt <- sample(x = unique(sort(d$time)), size = 10) 
d2 <- d
times <- sort(c(0,d$time))
seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)[c(1,2,5,12,90,125,200,241,267,268)]

test_that("conditional predictCSC with strata",{
  CSC.fitS <- CSC(Hist(time,event)~ strata(X1) + X5 + strata(X3) + X7 +X2,data=d, method = "breslow", surv.type = "survival")

  p1 <- predict(CSC.fitS, newdata = d2[1:10,], times = seqTime, landmark = 1, cause = 1, se = FALSE)
  #exportRes(p1$absRisk) 
  pGS <- rbind(c(NA, NA, NA, NA, 0.193715909198832, 0.258450268162697, 0.344406405530398, 0.405789224546312, NA, NA),
               c(NA, NA, NA, NA, 0.321903728009673, 0.42277637346952, 0.54718926951405, 0.628009479852877, NA, NA),
               c(NA, NA, NA, NA, 0.323271317692245, 0.423840829739195, 0.546860704948394, 0.62592540773196, NA, NA),
               c(NA, NA, NA, NA, 0.158112584306882, 0.212179101162182, 0.285911694623438, 0.340374382665936, NA, NA),
               c(NA, NA, NA, NA, 0.148123942798325, 0.19792013203262, 0.264500929712643, 0.312467653540823, NA, NA),
               c(NA, NA, NA, NA, 0.102115459721344, 0.174717894142052, 0.2985814768624, 0.2985814768624, NA, NA),
               c(NA, NA, NA, NA, 0.505837470210174, 0.586975852821309, 0.930987981351775, 0.930987981351775, NA, NA),
               c(NA, NA, NA, NA, 0.182897491614277, 0.242798381159347, 0.320509989781403, 0.37438654562716, NA, NA),
               c(NA, NA, NA, NA, 0.12077796085695, 0.207396556962911, 0.357254019375158, 0.357254019375158, NA, NA),
               c(NA, NA, NA, NA, 0.0798409036011334, 0.139105253559894, 0.247553338801094, 0.247553338801094, NA, NA))
 
  expect_equal(as.double(p1$absRisk),as.double(pGS), tolerance = 1e-8)
  expect_error(predict(CSC.fitS, newdata = d2[1:10,], times = seqTime, landmark = 1, cause = 1, se = TRUE))
})

# }}}

# }}}

# {{{ 5- Fast SE
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")
setkey(d,time)

m.CSC <- CSC(Hist(time, event) ~ X1+X6, data = d)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

test_that("iid minimal - no strata", {
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = TRUE, cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = TRUE, cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$absRisk.se,res3$absRisk.se)
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "minimal", average.iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$absRisk.se,res3$absRisk.se)
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    expect_equal(res2$absRisk.average.iid, t(apply(res3$absRisk.iid,2:3,mean)))
})

m.CSC <- CSC(Hist(time, event) ~ strata(X1)+X6, data = d)
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

test_that("iid minimal - strata", {
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = TRUE, cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = TRUE, cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$absRisk.se,res3$absRisk.se)    
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "minimal", average.iid = TRUE)
    res3 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    log.transform = FALSE, cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$absRisk.se,res3$absRisk.se)
    expect_equal(res1$absRisk.iid,res3$absRisk.iid)
    expect_equal(res2$absRisk.average.iid, t(apply(res3$absRisk.iid,2:3,mean)))
})


    
# }}}

# {{{ 6- Prediction when iid is stored
data(Melanoma, package = "riskRegression")
cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+strata(sex),
                          Hist(time,status)~age+strata(sex)),
             data=Melanoma)
cfit1$iid <- lapply(cfit1$model,iidCox)


res <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
               times=4,se=TRUE,band=TRUE)

# }}}

#----------------------------------------------------------------------
### test-predictCSC_vs_mstate.R ends here
