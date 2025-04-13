library(testthat)
context("Cause-specific Cox regression")
library(riskRegression)
library(rms)
library(survival)
library(prodlim)
library(data.table)
library(mstate)
library(prodlim)
library(timereg)
nsim.band <- 100
tmat <- trans.comprisk(2, names = c("0", "1", "2"))
## ** Data
## simulatenous events
df1 <- data.frame(time = rep(1:10,2),event = c(rep(1,10),rep(2,10)),X1 = 0:1)
df1$event1 <- as.numeric(df1$event == 1)
df1$event2 <- as.numeric(df1$event == 2)
set.seed(11)
dfS <- rbind(cbind(df1, grp = 1, X2 = rbinom(20,1,.4)),cbind(rbind(df1,df1), grp = 2, X2 = rbinom(20,1,.4)),cbind(df1, grp = 3, X2 = rbinom(20,1,.4)))
## distinct events
df2 <- data.frame(time = c(1:10,-0.1+1:10),event = c(rep(1,10),rep(2,10)),X1 = 0:1)
df2$event1 <- as.numeric(df2$event == 1)
df2$event2 <- as.numeric(df2$event == 2)

test_that("CSC vs prodlim",{
    data(Melanoma)
    nd <- data.frame(sex=factor(levels(Melanoma$sex)))
    A <- prodlim(Hist(time,status)~1,data=Melanoma)
    B <- CSC(Hist(time,status)~1,data=Melanoma)
    pB <- predictRisk(B,times=sort(unique(Melanoma$time)),newdata=nd[1,,drop=FALSE],cause=1)
    pA <- predictRisk(A,times=sort(unique(Melanoma$time)),newdata=nd[1,,drop=FALSE],cause=1)
    expect_equal(ignore_attr=TRUE,c(pB),pA,tolerance=1e3)
    a <- prodlim(Hist(time,status)~sex,data=Melanoma)
    b <- CSC(Hist(time,status)~strat(sex),data=Melanoma,fitter="cph")
    c <- CSC(Hist(time,status)~strat(sex),data=Melanoma,surv.type="survival",fitter="cph")
    pa <- predictRisk(a,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    pb <- predictRisk(b,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    pc <- predictRisk(c,times=c(0,10,100,1000,2000),newdata=nd,cause=1)
    expect_equal(ignore_attr=TRUE,c(pb),c(pa),tolerance=0.1)
    expect_equal(ignore_attr=TRUE,c(pb),c(pc),tolerance=0.00000001)
    u <- CSC(Hist(time,status)~strat(sex)+age+invasion,data=Melanoma,fitter="cph")
    v <- CSC(Hist(time,status)~strat(sex)+age+invasion,data=Melanoma,surv.type="survival",fitter="cph")
    pu <- predictRisk(u,times=c(0,10,100,1000,2000),newdata=Melanoma[c(17,84),],cause=1)
    pv <- predictRisk(v,times=c(0,10,100,1000,2000),newdata=Melanoma[c(17,84),],cause=1)
    expect_equal(ignore_attr=TRUE,c(pu),c(pv),tolerance=0.1)
    ## plot(a,cause=1,ylim=c(0,.5),confint=FALSE)
    ## lines(sort(unique(Melanoma$time)),pb[1,],lty=3)
    ## lines(sort(unique(Melanoma$time)),pb[2,],lty=3,col=2)
})

if (requireNamespace("pec",quietly=TRUE)){
    library(pec)
    test_that("predictSurv",{
        set.seed(17)
        d <- prodlim::SimSurv(100)
        f <- coxph(Surv(time,status)~X1+X2,data=d,x=TRUE,y=TRUE)
        h <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
        af <- predictRisk(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        bf <- 1-predictSurvProb(f,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        expect_equal(ignore_attr=TRUE,unname(af),unname(bf),tolerance = 1e-8)
        ah <- predictRisk(h,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        bh <- 1-predictSurvProb(h,newdata=d[c(17,88,3),],times=c(0,1,8.423,100,1000))
        colnames(bh) <- NULL
        expect_equal(ignore_attr=TRUE,unname(ah),unname(bh),tolerance = 1e-8)
    })
}
test_that("Cox models",{
    set.seed(17)
    d <- sampleData(100)
    a <- CSC(Hist(time,event)~X1+X2,data=d)
    A <- CSC(Hist(time,event)~X1+X2,data=d,surv.type="surv")
    a1 <- coxph(Surv(time,event==1)~X1+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~X1+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~X1+X2,data=d)
    expect_equal(ignore_attr=TRUE,coef(a$models[[1]]),coef(a1),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(a$models[[2]]),coef(a2),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(A$models[[2]]),coef(A2),tolerance = 1e-8)
})

test_that("strata",{
    set.seed(17)
    d <- sampleData(100)
    a <- CSC(Hist(time,event)~strata(X1)+X2,data=d)
    A <- CSC(Hist(time,event)~strata(X1)+X2,data=d,surv.type="surv")
    a1 <- coxph(Surv(time,event==1)~strata(X1)+X2,data=d)
    a2 <- coxph(Surv(time,event==2)~strata(X1)+X2,data=d)
    A2 <- coxph(Surv(time,event!=0)~strata(X1)+X2,data=d)
    expect_equal(ignore_attr=TRUE,coef(a$models[[1]]),coef(a1),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(a$models[[2]]),coef(a2),tolerance = 1e-8)
    expect_equal(ignore_attr=TRUE,coef(A$models[[2]]),coef(A2),tolerance = 1e-8)
})

test_that("strat and strata",{
    data(Melanoma)
    a <- CSC(Hist(time,status)~strat(sex)+age+invasion+logthick+strat(epicel)+strat(ulcer),data=Melanoma,fitter="cph")
    predictRisk(a,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    b <- CSC(Hist(time,status)~strata(sex)+age+invasion+logthick+strata(epicel)+strata(ulcer),data=Melanoma,fitter="coxph")
    pa <- predictRisk(a,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    pb <- predictRisk(b,times=c(0,100,1000,4000),newdata=Melanoma[c(17,77,188),],cause=2)
    expect_equal(ignore_attr=TRUE,pa,pb,tolerance=1e-6)
})





## ** No risk factor no strata
test_that("[predictCSC] (no strata, data=df1): compare to manual estimation",{
    CSC.exp <- CSC(Hist(time,event) ~ 1, data = df1)
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = FALSE)  
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE) 
    CSC.prodlim <- CSC(Hist(time,event) ~ 1, data = df1, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    expect_equal(ignore_attr=TRUE,pred.prodlim,pred.exp2)
    CSC.prodlimE <- CSC(Hist(time,event) ~ 1, data = df1, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  

    ## test baseline hazard
    lambda1 <- lambda2 <- 1/seq(20,2,by = -2)
    lambda2E <- predictCox(CSC.prodlimE$models[[2]], type = "hazard")$hazard
    expect_equal(ignore_attr=TRUE,predictCox(CSC.exp$models[[1]], type = "hazard")$hazard,
                 lambda1)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.exp$models[[2]], type = "hazard")$hazard,
                 lambda2)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.prodlim$models[[1]], type = "hazard")$hazard,
                 lambda1)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.prodlim$models[[2]], type = "hazard")$hazard,
                 lambda1+lambda2)
    expect_equal(ignore_attr=TRUE,basehaz(CSC.prodlimE$models[[2]])$hazard,
                 cumsum(lambda2E)
                 )

    ## test absolute risk
    survival <- c(1,exp(cumsum(-lambda1-lambda2)))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk),
                 cumsum(lambda1*survival))
    survival <- c(1,cumprod(1-(lambda1+lambda2)))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk),
                 cumsum(lambda1*survival))
    survival <- c(1,cumprod(1-lambda2E))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk),
                 cumsum(lambda1*survival))
})

test_that("[predictCSC] (no strata, data=df1): compare to manual estimation",{
    CSC.exp <- CSC(Hist(time,event) ~ 1, data = df2)
    pred.exp <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = FALSE)  
    pred.exp2 <- predict(CSC.exp, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    CSC.prodlim <- CSC(Hist(time,event) ~ 1, data = df2, surv.type = "survival", method = "breslow")
    pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    expect_equal(ignore_attr=TRUE,pred.exp2, pred.prodlim)
    CSC.prodlimE <- CSC(Hist(time,event) ~ 1, data = df2, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(NA), cause = 1, product.limit = TRUE)  
    ## test baseline hazard
    lambda1 <- 1/seq(19,1,by = -2)
    lambda2 <- 1/seq(20,2,by = -2)
    lambda2E <- as.data.table(predictCox(CSC.prodlimE$models[[2]], type = "hazard"))
    expect_equal(ignore_attr=TRUE,predictCox(CSC.exp$models[[1]], times = 1:10, type = "hazard")$hazard,
                 lambda1)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.exp$models[[2]], times = -0.1 + 1:10, type = "hazard")$hazard,
                 lambda2)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.prodlim$models[[1]], times = 1:10, type = "hazard")$hazard,
                 lambda1)
    expect_equal(ignore_attr=TRUE,predictCox(CSC.prodlim$models[[2]], type = "hazard")$hazard,
                 as.double(rbind(lambda2,lambda1)))
    expect_equal(ignore_attr=TRUE,basehaz(CSC.prodlimE$models[[2]])$hazard,
                 cumsum(lambda2E$hazard)
                 )
    ## test absolute risk
    vec.lambda <- cumsum(as.double(rbind(lambda2,lambda1)))
    survival <- exp(-matrix(vec.lambda, ncol = 2, byrow = TRUE)[,1])
    expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk),
                 cumsum(lambda1*survival))
    vec.lambda <- as.double(rbind(lambda2,lambda1))
    survival <- matrix(cumprod(1-vec.lambda), ncol = 2, byrow = TRUE)[,1]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk),
                 cumsum(lambda1*survival))
    vec.lambda <-lambda2E$hazard
    survival <- matrix(cumprod(1-vec.lambda), ncol = 2, byrow = TRUE)[,1]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk),
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
        expect_equal(ignore_attr=TRUE,pred.probtrans[,"pstate2"],
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
        expect_equal(ignore_attr=TRUE,surv.prodlim+as.double(haz1.prodlim$absRisk)+as.double(haz2.prodlim$absRisk),rep(1,nTimes))
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
    expect_equal(ignore_attr=TRUE,pred.exp2,pred.prodlim)
    
    CSC.prodlimE <- CSC(Hist(time,event) ~ strata(grp), data = dfS, surv.type = "survival", method = "efron")
    pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = data.frame(grp = 1:3), cause = 1, product.limit = TRUE)  

    ## baseline
    baseline1 <- as.data.table(predictCox(CSC.exp$models[[1]], type = "hazard"))
    lambda1 <- baseline1[strata=="grp=1",hazard]
    expect_equal(ignore_attr=TRUE,lambda1, 1/seq(20,2,by = -2))

    baseline2 <- as.data.table(predictCox(CSC.exp$models[[1]], type = "hazard"))
    lambda2 <- baseline2[strata=="grp=1",hazard]
    expect_equal(ignore_attr=TRUE,lambda2, 1/seq(20,2,by = -2))
        
    baseline1E <- as.data.table(predictCox(CSC.prodlimE$models[[1]], type = "hazard"))
    lambda1E.1 <- baseline1E[strata=="grp=1",hazard]
    lambda1E.2 <- baseline1E[strata=="grp=2",hazard]
    expect_equal(ignore_attr=TRUE,lambda1E.1, 1/seq(20,2,by = -2))

    baseline2E <- as.data.table(predictCox(CSC.prodlimE$models[[2]], type = "hazard"))
    lambda2E.1 <- baseline2E[strata=="grp=1",hazard]
    lambda2E.2 <- baseline2E[strata=="grp=2",hazard]
    ## 
        
    ## test absolute risk
    survival <- c(1,exp(cumsum(-lambda1-lambda2)))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk[1,]),
                 cumsum(lambda1*survival))
    expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk[1,]),
                 as.double(pred.exp$absRisk[2,]))
    expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk[1,]),
                 as.double(pred.exp$absRisk[3,]))

    survival <- c(1,cumprod(1-(lambda1+lambda2)))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk[1,]),
                 cumsum(lambda1*survival))
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk[1,]),
                 as.double(pred.prodlim$absRisk[2,]))
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk[1,]),
                 as.double(pred.prodlim$absRisk[3,]))

    survival <- c(1,cumprod(1-lambda2E.1))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=1"]),
                 cumsum(lambda1E.1*survival))
    survival <- c(1,cumprod(1-lambda2E.2))[1:10]
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk[pred.prodlimE$strata=="grp=2"]),
                 cumsum(lambda1E.2*survival))    
    expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk[1,]),
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
        expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk),
                     calcCIF(CSC.exp, X = as.matrix(newdata))
                     )
        pred.prodlim <- predict(CSC.prodlim, times = 1:10, newdata = newdata, cause = 1, product.limit = TRUE)
        expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk),
                     calcCIF(CSC.prodlim, X = as.matrix(newdata))
                     )
        pred.prodlimE <- predict(CSC.prodlimE, times = 1:10, newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk),
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
        expect_equal(ignore_attr=TRUE,as.double(pred.exp$absRisk),
                     calcCIF(CSC.exp, X = as.matrix(newdata))
                     )
    
        pred.prodlim <- predict(CSC.prodlim, times = sort(unique(df2$time)), newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(ignore_attr=TRUE,as.double(pred.prodlim$absRisk),
                     calcCIF(CSC.prodlim, X = as.matrix(newdata))
                     )

        pred.prodlimE <- predict(CSC.prodlimE, times = sort(unique(df2$time)), newdata = newdata, cause = 1, product.limit = TRUE)  
        expect_equal(ignore_attr=TRUE,as.double(pred.prodlimE$absRisk),
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
            expect_equal(ignore_attr=TRUE,pred.probtrans[,"pstate2"],
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

            expect_equal(ignore_attr=TRUE,surv.prodlim+as.double(haz1.prodlim$absRisk)+as.double(haz2.prodlim$absRisk),rep(1,nTimes))
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
             data = as.data.frame(d), keep = c("X1","X2","X16","Xcat2"),
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


    expect_equal(ignore_attr=TRUE,as.double(pred.RR1a$absRisk),pred.probtrans[,"pstate2"])
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1c$absRisk),pred.probtrans[,"pstate3"])
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1d$absRisk),pred.probtrans[,"pstate3"], tol = 5e-2)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR2a$absRisk),pred.probtrans[,"pstate2"])
    expect_equal(ignore_attr=TRUE,as.double(pred.RR2b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)

    # quantile(as.double(pred.RR1a$absRisk.se) - pred.probtrans[,"se2"])
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 5e-3)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 5e-3)
    expect_equal(ignore_attr=TRUE,as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)    
    expect_equal(ignore_attr=TRUE,as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 5e-3)    
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
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1a$absRisk),pred.probtrans[,"pstate2"])
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1b$absRisk),pred.probtrans[,"pstate2"], tol = 5e-3)
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1c$absRisk),pred.probtrans[,"pstate3"])
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1d$absRisk),pred.probtrans[,"pstate3"], tol = 1e-2)
        expect_equal(ignore_attr=TRUE,as.double(pred.RR2a$absRisk),pred.probtrans[,"pstate2"], tol = 1e-2)
        expect_equal(ignore_attr=TRUE,as.double(pred.RR2b$absRisk),pred.probtrans[,"pstate2"], tol = 1e-1)
        #if(iX==0){
        # quantile(as.double(pred.RR1a$absRisk.se) - pred.probtrans[,"se2"])
        # quantile(as.double(pred.RR1c$absRisk.se) - pred.probtrans[,"se3"])
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        expect_equal(ignore_attr=TRUE,as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        if(iX==0){
            expect_equal(ignore_attr=TRUE,as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 1e-2)
            expect_equal(ignore_attr=TRUE,as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 1e-2)
        }
        expect_equal(ignore_attr=TRUE,as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)    
        expect_equal(ignore_attr=TRUE,as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-2)
        ## }else{
        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR1a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)

        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR1b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)
        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR1c$absRisk.se),pred.probtrans[,"se3"], tol = 1e-1)
        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR1d$absRisk.se),pred.probtrans[,"se3"], tol = 1e-1)
        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR2a$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)    
        ##     expect_equal(ignore_attr=TRUE,as.double(pred.RR2b$absRisk.se),pred.probtrans[,"se2"], tol = 1e-1)
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

    expect_equal(ignore_attr=TRUE,pred.RR1a_strata$absRisk,pred.RR1a$absRisk)
    expect_equal(ignore_attr=TRUE,pred.RR1b_strata$absRisk,pred.RR1b$absRisk)
    expect_equal(ignore_attr=TRUE,pred.RR1c_strata$absRisk,pred.RR1c$absRisk)
    expect_equal(ignore_attr=TRUE,pred.RR1d_strata$absRisk,pred.RR1d$absRisk)

    CSC.RR2_strata<- CSC(Hist(time,event)~X1+X2+X16+strata(grp), data = d2, surv.type = "survival", method = "breslow")
    pred.RR2a_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                         keep.newdata = FALSE, se = TRUE, product.limit = TRUE)
    pred.RR2b_strata <- predict(CSC.RR2_strata, cbind(newdata,grp=1), cause = 1, time = pred.probtrans[,"time"],
                                keep.newdata = FALSE, se = TRUE, product.limit = FALSE)

    expect_equal(ignore_attr=TRUE,pred.RR2a_strata$absRisk,pred.RR2a$absRisk)
    expect_equal(ignore_attr=TRUE,pred.RR2b_strata$absRisk,pred.RR2b$absRisk)
})


test_that("predict.CSC (Melanoma): compare to mstate",{
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
    ## mstate
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
    expect_equal(ignore_attr=TRUE,as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR1$absRisk))
    pred.RR2 <- predict(cfit1, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
    expect_equal(ignore_attr=TRUE,as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR2$absRisk), tol = 5e-3)
    pred.RR3 <- predict(cfit2, newdata = newdata, time = pred.probtrans[,"time"], cause = 1, product.limit = FALSE) 
    expect_equal(ignore_attr=TRUE,as.double(pred.probtrans[,"pstate2"]),as.double(pred.RR3$absRisk), tol = 1e-1)
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
  expect_equal(ignore_attr=TRUE,predC_auto$absRisk,predC_manuel)
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
  expect_equal(ignore_attr=TRUE,predC_auto$absRisk,predC_manuel)
  # predC_auto$absRisk-predC_manuel
})

test_that("[predictCSC] diag no strata", {
    ## * [predictCSC] Argument diag
    cat("[predictCSC] Argument \'diag\' \n")
    set.seed(10)
    dt <- sampleData(75, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]
    e.CSC <- CSC(Hist(time, event) ~ X1*X6, data = dt)
    GS <- predict(e.CSC, newdata = dt, times = dt$time, se = TRUE, iid = TRUE, average.iid = TRUE, cause = 1)
    test <- predict(e.CSC, newdata = dt, times = dt$time, se = TRUE, iid = TRUE, average.iid = TRUE, diag = TRUE, cause = 1)
    test2 <- predict(e.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE, cause = 1)
    ## estimates
    expect_equal(ignore_attr=TRUE,dt$time, as.double(test$time))
    expect_equal(ignore_attr=TRUE,diag(GS$absRisk), as.double(test$absRisk))
    ## se
    expect_equal(ignore_attr=TRUE,diag(GS$absRisk.se), test$absRisk.se[,1])
    ## iid
    GS.iid.diag <- do.call(cbind,lapply(1:NROW(dt),
                                        function(iN){GS$absRisk.iid[,iN,iN]}))
    expect_equal(ignore_attr=TRUE,GS.iid.diag, test$absRisk.iid[,1,])
    ## average.iid
    expect_equal(ignore_attr=TRUE,rowMeans(GS.iid.diag), test$absRisk.average.iid[,1])
    expect_equal(ignore_attr=TRUE,test$absRisk.average.iid, test2$absRisk.average.iid)
    ## average.iid with factor - diag=FALSE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test3 <- predict(e.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE, cause = 1)
    expect_equal(ignore_attr=TRUE,5*GS$absRisk.average.iid, test3$absRisk.average.iid[[1]])
    expect_equal(ignore_attr=TRUE,apply(GS$absRisk.iid, 1:2, function(x){sum(x * (1:length(dt$time)))/length(x)}),
                 test3$absRisk.average.iid[[2]])
    ## average.iid with factor - diag=FALSE, time varying factor
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(rnorm(NROW(dt)*length(dt$time)), nrow = NROW(dt), ncol = length(dt$time)))
    test4 <- predict(e.CSC, newdata = dt, times = dt$time, cause = 1,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)
    ## iObs <- 1
    expect_equal(ignore_attr=TRUE,do.call(rbind,lapply(1:NROW(dt), function(iObs){rowMeans(GS$absRisk.iid[iObs,,] * t(attr(average.iid,"factor")[[1]]))})),
                 test4$absRisk.average.iid[[1]])
    ## average.iid with factor - diag=TRUE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1, byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test5 <- predict(e.CSC, newdata = dt, times = dt$time, cause = 1,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = TRUE)
    expect_equal(ignore_attr=TRUE,5*test$absRisk.average.iid, test5$absRisk.average.iid[[1]])
    expect_equal(ignore_attr=TRUE,rowMeans(rowMultiply_cpp(GS.iid.diag, 1:length(dt$time))),
                 test5$absRisk.average.iid[[2]][,1])
})


test_that("[predictCSC] diag strata", {
    set.seed(10)
    dt <- sampleData(75, outcome = "competing.risks")[,.(time,event,X1,X2,X6)]
    eS.CSC <- CSC(Hist(time, event) ~ strata(X1) + X6, data = dt)
    GS <- predict(eS.CSC, newdata = dt, times = dt$time, se = TRUE, iid = TRUE, average.iid = TRUE, cause = 1)
    test <- predict(eS.CSC, newdata = dt, times = dt$time,
                    se = TRUE, iid = TRUE, average.iid = TRUE, diag = TRUE, cause = 1)
    test2 <- predict(eS.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = TRUE, diag = TRUE, cause = 1)
    ## estimates
    expect_equal(ignore_attr=TRUE,dt$time, as.double(test$time))
    expect_equal(ignore_attr=TRUE,diag(GS$absRisk), as.double(test$absRisk))
    ## se
    expect_equal(ignore_attr=TRUE,diag(GS$absRisk.se), test$absRisk.se[,1])
    ## iid
    GS.iid.diag <- do.call(cbind,lapply(1:NROW(dt),
                                        function(iN){GS$absRisk.iid[,iN,iN]}))
    expect_equal(ignore_attr=TRUE,GS.iid.diag, test$absRisk.iid[,1,])
    ## average.iid
    expect_equal(ignore_attr=TRUE,rowMeans(GS.iid.diag), test$absRisk.average.iid[,1])
    expect_equal(ignore_attr=TRUE,test$absRisk.average.iid, test2$absRisk.average.iid)
    ## average.iid with factor - diag=FALSE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test3 <- predict(eS.CSC, newdata = dt, times = dt$time,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE, cause = 1)
    expect_equal(ignore_attr=TRUE,5*GS$absRisk.average.iid, test3$absRisk.average.iid[[1]])
    expect_equal(ignore_attr=TRUE,apply(GS$absRisk.iid, 1:2, function(x){sum(x * (1:length(dt$time)))/length(x)}),
                 test3$absRisk.average.iid[[2]])
    ## average.iid with factor - diag=FALSE, time varying factor
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(rnorm(NROW(dt)*length(dt$time)), nrow = NROW(dt), ncol = length(dt$time)))
    test4 <- predict(eS.CSC, newdata = dt, times = dt$time, cause = 1,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = FALSE)
    expect_equal(ignore_attr=TRUE,do.call(rbind,lapply(1:NROW(dt), function(iObs){rowMeans(GS$absRisk.iid[iObs,,] * t(attr(average.iid,"factor")[[1]]))})),
                 test4$absRisk.average.iid[[1]])
    ## average.iid with factor - diag=TRUE
    average.iid <- TRUE
    attr(average.iid,"factor") <- list(matrix(5, nrow = NROW(dt), ncol = 1, byrow = TRUE),
                                       matrix(1:NROW(dt), nrow = NROW(dt), ncol = 1))
    test5 <- predict(eS.CSC, newdata = dt, times = dt$time, cause = 1,
                     se = FALSE, iid = FALSE, average.iid = average.iid, diag = TRUE)
    expect_equal(ignore_attr=TRUE,5*test$absRisk.average.iid, test5$absRisk.average.iid[[1]])
    expect_equal(ignore_attr=TRUE,rowMeans(rowMultiply_cpp(GS.iid.diag, 1:length(dt$time))),
                 test5$absRisk.average.iid[[2]][,1])
})


test_that("iid average - non and semi parametric (hazard)", {
    ## * [predictCSC] Average iid
    cat("[predictCSC] Average iid \n")
    set.seed(10)
    d <- sampleData(100, outcome = "competing.risks")
    obs.firstEvent <- which.min(d$time)
    strata.firstEvent <- 3
    seqTime <- unique(c(0,min(d$time),d$time[sample.int(n = 100, size = 10)],1e8))
    index.firstEvent <- 2
    ## ** non parametric
    for(iType in c("hazard","survival")){ ## iType <- "hazard"
        m.CSC <- CSC(Hist(time, event) ~ strata(X1) + strata(X2),
                     data = d, surv.type = iType)
        res1 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, iid = TRUE, average.iid = TRUE)
        res2 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, average.iid = TRUE)
        expect_equal(ignore_attr=TRUE,res1$absRisk.average.iid,res2$absRisk.average.iid)
        expect_true(all(res1$absRisk.average.iid[,1]==0))
        expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))
        expect_equal(ignore_attr=TRUE,apply(res1$absRisk.iid,1:2,mean),
                     res2$absRisk.average.iid)
        ## compare to fixed value    
        ## d[time==min(time),]
        ## levels(predictCox(m.CSC$models[[1]])$strata)
        expect_equal(ignore_attr=TRUE,res1$absRisk.iid[, index.firstEvent,obs.firstEvent],
                     iidCox(m.CSC$models[[1]], return.object = FALSE)$IFhazard[[strata.firstEvent]][,1])
    }
    ## ** semi parametric
    ## GIVES WARNING
    ## Estimated risk outside the range [0,1].
    for(iType in c("hazard","survival")){ ## iType <- "hazard"
        m.CSC <- CSC(Hist(time, event) ~ X1*X6 + strata(X2), data = d, surv.type = iType)
        suppressWarnings(res1 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, iid = TRUE, average.iid = TRUE))
        suppressWarnings(res2 <- predict(m.CSC, times = seqTime, newdata = d,
                        cause = 1, average.iid = TRUE))
        GS3 <- res1$absRisk.average.iid
        GS23 <- res1$absRisk.average.iid
        expect_equal(ignore_attr=TRUE,res1$absRisk.average.iid,res2$absRisk.average.iid)
        expect_true(all(res1$absRisk.average.iid[,1]==0))
        expect_true(all(is.na(res1$absRisk.average.iid[,length(seqTime)])))
        expect_equal(ignore_attr=TRUE,apply(res1$absRisk.iid,1:2,mean),
                     res2$absRisk.average.iid)
    }
})


test_that("predictSurv (type=survival)", {
    ## * [predictCSC] survival
    cat("[predictCSC] survival \n")
    ## ** Simulate data
    set.seed(10)
    d <- sampleData(2e2, outcome = "competing.risk")
    d.pred <- d[5:15]
    seqTime <- c(0,d[["time"]][5:15],2.45,1e5)
    ## ** type=survival CSC
    e.CSC <- CSC(Hist(time, event)~ strata(X1) + strata(X2), data = d, surv.type = "hazard")
    e.cox <- coxph(Surv(time, event>0)~ strata(X1) + strata(X2), data = d, x = TRUE, y = TRUE)
    ## e.cox <- coxph(Surv(time, event>0)~ 1, data = d, x = TRUE, y = TRUE)
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime,
                    product.limit = FALSE, iid = TRUE, se = TRUE, average.iid = TRUE,
                    store.iid = "minimal", keep.newdata = FALSE)
    GS <- predictCox(e.cox, newdata = d.pred, times = seqTime, type = "survival",
                     iid = TRUE, se = TRUE, average.iid = TRUE)
    expect_equal(ignore_attr=TRUE,GS$survival, test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    testPL <- predict(e.CSC, type = "survival", newdata = d.pred, times = seqTime,
                      product.limit = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE,
                      store.iid = "minimal")
    GSPL <- predictCoxPL(e.cox, newdata = d.pred, times = seqTime,
                         iid = TRUE, se = TRUE, average.iid = TRUE)
    ratioSurv <- test$survival/testPL$survival
    testPL$survival.iid2 <- testPL$survival.iid
    for(iObs in 1:NROW(d)){
        testPL$survival.iid2[iObs,,] <- testPL$survival.iid[iObs,,]*t(ratioSurv)
    }
    expect_equal(ignore_attr=TRUE,GSPL$survival,testPL$survival)
    expect_equal(ignore_attr=TRUE,GSPL$survival.se[!is.infinite(ratioSurv)], (testPL$survival.se * ratioSurv)[!is.infinite(ratioSurv)])
    expect_equal(ignore_attr=TRUE,GSPL$survival.iid[!is.nan(testPL$survival.iid2)],testPL$survival.iid2[!is.nan(testPL$survival.iid2)])
    expect_equal(ignore_attr=TRUE,apply(testPL$survival.iid,1:2,mean),
                 testPL$survival.average.iid)
    # "predictSurv (type=survival,diag)"
    test <- predict(e.CSC, type = "survival", newdata = d.pred, times = d.pred$time,
                    diag = TRUE, product.limit = FALSE, iid = TRUE, se = TRUE, average.iid = TRUE,
                    store.iid = "minimal")
    GS <- predictCox(e.cox, newdata = d.pred, times = d.pred$time,
                     diag = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE)
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    testPL <- predict(e.CSC, type = "survival", newdata = d.pred, times = d.pred$time,
                      diag = TRUE, product.limit = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE,
                      store.iid = "minimal")
    suppressWarnings(GSPL <- predictCoxPL(e.cox, newdata = d.pred, times = d.pred$time,
                         diag = TRUE, iid = TRUE, se = TRUE, average.iid = TRUE))
    ratioSurv <- test$survival/testPL$survival
    testPL$survival.iid2 <- testPL$survival.iid
    for(iObs in 1:NROW(d)){
        testPL$survival.iid2[iObs,,] <- testPL$survival.iid[iObs,,]*t(ratioSurv)
    }
    expect_equal(ignore_attr=TRUE,GSPL$survival,testPL$survival)
    expect_equal(ignore_attr=TRUE,GSPL$survival.se[!is.infinite(ratioSurv)], (testPL$survival.se * ratioSurv)[!is.infinite(ratioSurv)])
    expect_equal(ignore_attr=TRUE,GSPL$survival.iid[!is.nan(testPL$survival.iid2)],testPL$survival.iid2[!is.nan(testPL$survival.iid2)])
    expect_equal(ignore_attr=TRUE,apply(testPL$survival.iid,1:2,mean),
                 testPL$survival.average.iid)
  # "[predictCSC] vs. predictCox (no strata) - surv.type=\"survival\"",{
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE, cause = 1)
    GS <- predictCox(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE))
    suppressWarnings(GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
})
test_that("[predictCSC] vs. predictCox (no strata) - surv.type=\"survival\"",{
    set.seed(8)
    d = sampleData(74)
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE, cause = 1))
    suppressWarnings(GS <- predictCox(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE))
    suppressWarnings(GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]],times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    # test_that("[predictCSC] vs. predictCox (strata) - surv.type=\"survival\"",{
    e.CSC <- CSC(Hist(time, event)~ X6 + strata(X1), data = d, surv.type = "survival")
    jumpTime <- e.CSC$eventTimes[e.CSC$eventTimes <= max(seqTime)]
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE))
    suppressWarnings(GS <- predictCox(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE))
    suppressWarnings(GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    ## different strata for each cause
    e.CSC <- CSC(list(Hist(time, event)~ X6 + strata(X1),
                      Hist(time, event)~ X6 + strata(X2)),
                 data = d, surv.type = "survival")
    test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = FALSE, iid = TRUE, average.iid = TRUE, se = TRUE)
    GS <- predictCox(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                     newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE)
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    suppressWarnings(test <- predict(e.CSC, type = "survival", times = jumpTime,
                    newdata = d.pred, product.limit = TRUE, iid = TRUE, average.iid = TRUE, se = TRUE))
    suppressWarnings(GS <- predictCoxPL(e.CSC$models[["OverallSurvival"]], type = "survival", times = jumpTime,
                       newdata = d.pred, iid = TRUE, average.iid = TRUE, se = TRUE))
    expect_equal(ignore_attr=TRUE,GS$survival,test$survival)
    expect_equal(ignore_attr=TRUE,GS$survival.se,test$survival.se)
    expect_equal(ignore_attr=TRUE,GS$survival.iid,test$survival.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
})

test_that("[predictCSC] for survival surv.type=\"survival\" (internal consistency)",{
    ## ** type=hazard
    tau <- median(d$time)
    e.CSC <- CSC(Hist(time, event)~ X6, data = d, surv.type = "hazard")
    GS <- predict(e.CSC, type = "survival", times = tau,
                  newdata = d, product.limit = FALSE, iid = TRUE, average.iid = TRUE,
                  store.iid = "minimal")
    test <- predict(e.CSC, type = "survival", times = tau,
                    newdata = d, product.limit = FALSE, iid = FALSE, average.iid = TRUE,
                    store.iid = "minimal")
    factor <- TRUE
    attr(factor,"factor") <- list(matrix(5, nrow = NROW(d), ncol = 1),
                                  matrix(1:NROW(d), nrow = NROW(d), ncol = 1))
    test2 <- predict(e.CSC, type = "survival", times = tau, se = FALSE,
                     newdata = d, product.limit = FALSE, iid = FALSE, average.iid = factor,
                     store.iid = "minimal")
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid,test$survival.average.iid)
    expect_equal(ignore_attr=TRUE,GS$survival.average.iid[,1],test$survival.average.iid[,1])
    expect_equal(ignore_attr=TRUE,5*GS$survival.average.iid[,1],test2$survival.average.iid[[1]][,1])
    expect_equal(ignore_attr=TRUE,rowMeans(rowMultiply_cpp(GS$survival.iid[,1,],1:NROW(d))),
                 test2$survival.average.iid[[2]][,1])
})
test_that("Prediction with CSC - sorted vs. unsorted times",{
    ## ** Order of the prediction times
    cat("[predictCSC] order of the prediction times \n")
    data(Melanoma, package = "riskRegression")
    times2 <- sort(c(0,0.9*min(Melanoma$time),Melanoma$time[5],max(Melanoma$time)*1.1))
    newOrder <- sample.int(length(times2),length(times2),replace = FALSE)
    fit.CSC <- CSC(Hist(time,status) ~ thick, data = Melanoma)
    predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1, keep.times = FALSE)
    predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1, keep.times = FALSE)
    expect_equal(ignore_attr=TRUE,predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
    # Prediction with CSC (strata) - sorted vs. unsorted times"
    fit.CSC <- CSC(Hist(time,status) ~ thick + strat(invasion), data = Melanoma)
    predictionUNS <- predict(fit.CSC, times = times2[newOrder], newdata = Melanoma, cause = 1)
    predictionS <- predict(fit.CSC, times = times2, newdata = Melanoma, cause = 1)
    expect_equal(ignore_attr=TRUE,predictionS$absRisk, predictionUNS$absRisk[,order(newOrder)])
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
    expect_equal(ignore_attr=TRUE,out1$survival.average.iid,out2$survival.average.iid, tol = 1e-8)
})

## ** Prediction when iid is stored
cat("[predictCSC] prediction after storing iid \n")
test_that("[predictCSC] prediction after storing iid", {
    data(Melanoma, package = "riskRegression")
    cfit1 <- CSC(formula=list(Hist(time,status)~age+logthick+epicel+strata(sex),
                              Hist(time,status)~age+strata(sex)),
                 data=Melanoma)
    GS <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
                  times=4,se=TRUE,band=TRUE)
    cfit1 <- iidCox(cfit1)
    res <- predict(cfit1,newdata=Melanoma[1,,drop=FALSE],cause=1,
                   times=4,se=TRUE,band=TRUE)
    expect_equal(ignore_attr=TRUE,GS,res)
})
test_that("[predictCSC] (no strata): NA after last event",{
    ## *** no strata
    data(Melanoma, package = "riskRegression")
    times1 <- unique(Melanoma$time)
    times2 <- c(0,0.9*min(times1),times1*1.1)
    dataset1 <- Melanoma[sample.int(n = nrow(Melanoma), size = 12),]
    fit.CSC <- CSC(Hist(time,status) ~ thick*age, data = Melanoma, fitter = "cph")
    test.times <- max(Melanoma$time) + c(-1e-1,0,1e-1)
    prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
    expect_equal(ignore_attr=TRUE,as.vector(is.na(prediction$absRisk)), c(FALSE, FALSE, TRUE))
    # ("[predictCSC] (no strata): no event before prediction time",{
    test.times <- min(Melanoma$time)-1e-5
    prediction <- predict(fit.CSC, times = test.times, newdata = Melanoma[1,,drop = FALSE], cause = 1)
    expect_equal(ignore_attr=TRUE,as.double(prediction$absRisk), 0)
    # ("[predictCSC] (no strata): beyond last event time / negative timepoints / NA in timepoints",{
    test.times.1 <- c(10,3000,5000)
    newd <- Melanoma
    newd <- Melanoma[1,,drop=FALSE]
    prediction.1 <- predict(fit.CSC, times = test.times.1, newdata = newd, cause = 1)
    test.times.2 <- c(10,300,5000)
    prediction.2 <- predict(fit.CSC, times = test.times.2, newdata = newd, cause = 1)
    expect_equal(ignore_attr=TRUE,prediction.1$absRisk[,3],prediction.2$absRisk[,3])
    expect_equal(ignore_attr=TRUE,unname(predict(fit.CSC, times = -1, newdata = newd, cause = 1)$absRisk),
                 matrix(0,nrow = nrow(newd), ncol = 1))
    expect_error(predict(fit.CSC, times = c(test.times,NA), newdata = newd, cause = 1))
})

test_that("[predictCSC](strata): NA after last event",{
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
    for(Ttempo in 1:nrow(dt.times)){
        test.times <- sort(unlist(dt.times[Ttempo, .(beforeLastTime, LastTime, afterLastTime)]))
    prediction <- predict(fit.CSC, times = test.times, newdata = data.test, cause = 1)
    expect_equal(ignore_attr=TRUE,unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(ignore_attr=TRUE,unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
    expect_equal(ignore_attr=TRUE,unname(is.na(prediction$absRisk[Ttempo,])), c(FALSE, FALSE, TRUE))
  }
})


test_that("After last event - fully stratified model",{
    set.seed(10)
    d <- sampleData(1e2)
    setkeyv(d, "time")
    d[X1==0,event := c(event[1:(.N-1)],0)]
    d[X1==1,event := c(event[1:(.N-1)],1)]
    tau <- c(d[,max(time),by="X1"][[2]],10000)
    ## one strata variable
    X <- unique(d[,"X1",drop=FALSE])
    e.CSC <- CSC(Hist(time,event) ~ strata(X1), data = d)
    ## plot(prodlim(Hist(time,event) ~ X1, data = d))
    test1 <- as.data.table(predict(e.CSC, times = tau, newdata = X, se = TRUE, cause = 1))
    test2 <- as.data.table(predict(e.CSC, times = tau, newdata = X, se = TRUE, cause = 2))
    expect_equal(ignore_attr=TRUE,test1[strata == "X1=1" & times == 10000,absRisk], test1[strata == "X1=1" & times == d[X1==1,max(time)],absRisk], tol = 1e-6)
    expect_equal(ignore_attr=TRUE,test1[strata == "X1=1" & times == 10000,absRisk] + test2[strata == "X1=1" & times == 10000,absRisk],1,tol = 1e-6)
    expect_equal(ignore_attr=TRUE,test1[strata == "X1=1" & times == 10000,absRisk.se], test1[strata == "X1=1" & times == d[X1==1,max(time)],absRisk.se], tol = 1e-6)
    expect_true(is.na(test1[strata == "X1=0" & times == 10000,absRisk]))
    expect_true(is.na(test1[strata == "X1=0" & times == 10000,absRisk.se]))
    ## two strata variables
    X <- unique(d[,c("X1","X2"),drop=FALSE])
    e2.CSC <- CSC(Hist(time,event) ~ strata(X1)+strata(X2), data = d)
    ## plot(prodlim(Hist(time,event) ~ X1 + X2, data = d),cause =2)
    test1 <- as.data.table(predict(e2.CSC, times = tau, newdata = X, se = TRUE, cause = 1))
    test2 <- as.data.table(predict(e2.CSC, times = tau, newdata = X, se = TRUE, cause = 2))
    expect_equal(ignore_attr=TRUE,test1[strata != "X1=0 X2=0" & times == 10000,absRisk] + test2[strata != "X1=0 X2=0" & times == 10000,absRisk], rep(1,3), tol = 1e-6)
    expect_true(all(is.na(test1[strata == "X1=0 X2=0" & times == 10000,absRisk])))
    expect_true(all(is.na(test1[strata == "X1=0 X2=0" & times == 10000,absRisk.se])))
    expect_true(all(is.na(test2[strata == "X1=0 X2=0" & times == 10000,absRisk])))
    expect_true(all(is.na(test2[strata == "X1=0 X2=0" & times == 10000,absRisk.se])))
})

test_that("[predictCSC]: iid minimal", {
    ## ** Argument store.iid
    cat("[predictCSC] Argument \'store.iid\' \n")
    set.seed(10)
    d <- sampleData(1e2, outcome = "competing.risks")
    setkey(d,time)
    m.CSC <- CSC(Hist(time, event) ~ X1+X6, data = d)
    seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
    newdata <- d
    ## head(newdata[,.(time,event,X1,X6)])
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE, average.iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(ignore_attr=TRUE,res1$absRisk.se,res2$absRisk.se)
    expect_equal(ignore_attr=TRUE,res1$absRisk.iid,res2$absRisk.iid)
    expect_equal(ignore_attr=TRUE,res1$absRisk.average.iid, apply(res2$absRisk.iid,1:2,mean))
    m2.CSC <- iidCox(m.CSC, store.iid = "minimal")
    res1bis <- predict(m2.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    se = TRUE, iid = TRUE, average.iid = TRUE)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.se,res2$absRisk.se)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.iid,res2$absRisk.iid)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.average.iid, apply(res2$absRisk.iid,1:2,mean))
    # with strata
    m.CSC <- CSC(Hist(time, event) ~ strata(X1)+X6, data = d)
    res1 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "minimal", se = TRUE, iid = TRUE, average.iid = TRUE)
    res2 <- predict(m.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(ignore_attr=TRUE,res1$absRisk.se,res2$absRisk.se)
    expect_equal(ignore_attr=TRUE,res1$absRisk.iid,res2$absRisk.iid)
    expect_equal(ignore_attr=TRUE,res1$absRisk.average.iid, apply(res2$absRisk.iid,1:2,mean))
    m2.CSC <- iidCox(m.CSC, store.iid = "minimal")
    res1bis <- predict(m2.CSC, times = seqTime, newdata = newdata,
                    cause = 1,
                    se = TRUE, iid = TRUE, average.iid = TRUE)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.se,res2$absRisk.se)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.iid,res2$absRisk.iid)
    expect_equal(ignore_attr=TRUE,res1bis$absRisk.average.iid, apply(res2$absRisk.iid,1:2,mean))
})


test_that("[predictCSC]: surv.type", {
    data(Melanoma)
    Melanoma$id=1:nrow(Melanoma)
    Melanoma$status[196:205]=3
    fit2 <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma,
                surv.type="surv",cause=2)
    r2 <- predictRisk(fit2,cause=2,newdata=Melanoma,times=1000)
    Melanoma$status[196:205]=1
    fitGS <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma,
                 surv.type="surv",cause=2)
    rGS <- predictRisk(fitGS,cause=2,newdata=Melanoma,times=1000)
    expect_equal(ignore_attr=TRUE,r2,rGS, tol = 1e-6)
})


## *** tentative fix
test_that("product.limit=-1",{
    ## ** from Vi Friday 22-07-29 at 13:09
    simData <- function(n, sd = 1){
        m <- lava::lvm(~ X)
        lava::distribution(m, ~ X) <- lava::normal.lvm(mu = 0, sd = sd)
        lava::transform(m, X_2 ~ X) <- function(x){x ^ 2}
        lava::distribution(m, ~ Y) <- lava::binomial.lvm(p = 0.25)
        lava::distribution(m, ~ t0) <- lava::coxWeibull.lvm(scale = 2/100)
        lava::distribution(m, ~ t1) <- lava::coxWeibull.lvm(scale = 1/100)
        lava::distribution(m, ~ t2) <- lava::coxWeibull.lvm(scale = 0.5/100)
        lava::regression(m, t1 ~ X + X_2 + Y + X * Y) <- c(0.75, 1, 0.5, -0.1)
        lava::regression(m, t2 ~ X) <- c(0.25)
        m <- lava::eventTime(m, time ~ min(t0 = 0, t1 = 1, t2 = 2), 'status')
        simdat <- lava::sim(m, n)
        data.table::setDT(simdat)
        simdat[, -c('t0', 't1', 't2', 'X_2')]
    }
    ## set.seed(1)
    ## dt <- simData(50, sd = 1)[,.(time,status,X,Y)]
    dt <- data.table("time" = c(0.614, 0.704, 0.71, 0.758, 0.936, 1.11, 1.13, 1.147, 1.384, 1.586, 1.727, 1.804, 1.865, 1.922, 2.026, 2.04, 2.165, 2.229, 2.444, 2.639, 2.725, 2.772, 2.833, 2.927, 2.928, 3.052, 3.067, 3.13, 3.183, 3.411, 3.693, 3.796, 4.098, 4.348, 4.491, 4.565, 4.766, 4.834, 4.88, 5.124, 5.555, 5.609, 6.009, 7.009, 7.183, 7.239, 7.741, 8.254, 8.418, 10.579), 
                     "status" = c(1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 1, 2, 1, 1, 0, 0, 0, 1, 1, 1, 0, 2, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 0), 
                     "X" = c(1.764, -1.905, 1.324, 1.474, -2.592, 1.578, 0.596, 1.111, 0.38, -1.919, 0.865, 0.678, 1.092, 0.677, 0.568, -0.193, -0.43, 0.278, -0.156, 1.314, 0.045, 0.763, -0.389, 1.593, -1.168, 0.345, 0.164, -0.923, 1.155, -0.11, -2.129, -1.363, -0.169, -0.715, -0.811, -0.069, -0.146, -0.573, -0.008, 0.612, -0.823, -0.636, -0.164, -0.057, 0.307, -1.174, 0.616, 0.129, -0.924, -0.195))
    dt$X2 <- dt$X^2
    dt[,c("X","X2") := .(X-mean(X),X2-mean(X2))]
    index.obs <- 5
    dtPred <- dt[index.obs]
    ## ggplot(dt,aes(time,X)) + geom_point()
    e.CSC <- CSC(Hist(time, status) ~ X + X2, data = dt)
    ## *** issue
    if(FALSE){
        ls.beta <- coef(e.CSC)
        e.riskPL <- predict(e.CSC, newdata = dtPred, times = 5, cause = 1, product.limit = TRUE)$absRisk
        e.riskPL
        e.survPL <- predictCoxPL(e.CSC$models[[1]], newdata = dtPred, times = 5)$survival
        e.survPL
        eXb1 <- exp(cbind(dt$X,dt$X2) %*% ls.beta[[1]])
        lambda01 <-  (dt$status==1) / rev(cumsum(rev(eXb1))) ## one hazard is bigger than 1!
        ## range(predictCox(e.CSC$models[[1]], type = "hazard")$hazard - lambda01)
        dt$status[49]/sum(eXb1[49:50])
        e.riskEXP <- predict(e.CSC, newdata = dtPred, times = dt$time, cause = 1, product.limit = FALSE)$absRisk
        e.riskEXP
        eXb2 <- exp(cbind(dt$X,dt$X2) %*% ls.beta[[2]])
        lambda02 <- (dt$status==2) / rev(cumsum(rev(eXb2)))
        ## range(predictCox(e.CSC$models[[2]], type = "hazard")$hazard - lambda02)
        St <- exp(-cumsum(eXb1[index.obs]*lambda01-eXb2[index.obs]*lambda02))
        ## range(St - predict(e.CSC, newdata = dtPred, type = "survival", product.limit = FALSE, times = dt$time)$survival)
        StM1 <- c(1,St[1:(length(St)-1)])
        risk <- cumsum(lambda01*eXb1[index.obs] * StM1)
        range(risk - e.predEXP$absRisk)
    }
    eFIX.riskEXP <- predict(e.CSC, newdata = dtPred, times = dt$time, cause = 1,
                            se = TRUE, iid = TRUE, product.limit = -1)
    GS <- structure(c(0.132050275998517, 0.274110503238186, 0.397566794786161, 
                      0.508203384443532, 0.610377243723005, 0.610377243723005, 0.610377243723005, 
                      0.73909177290503, 0.840418036467214, 0.840418036467214, 0.840418036467214, 
                      0.921770849343348, 0.921770849343348, 0.921770849343348, 0.986517721864252, 
                      0.986517721864252, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), dim = c(1L, 50L))
    ## cap risks at 1
    expect_equal(ignore_attr=TRUE,eFIX.riskEXP$absRisk, GS, tol = 1e-6)
    expect_true(max(eFIX.riskEXP$absRisk)<=1)
    ## corresponding se=0
    expect_true(all(abs(eFIX.riskEXP$absRisk.se[eFIX.riskEXP$absRisk==1])<1e-6))
    ## corresponding IF=0
    expect_true(all(abs(eFIX.riskEXP$absRisk.iid[,eFIX.riskEXP$absRisk==1,1])<1e-6))
    ## check correct  averaging (via C++)
    eALL.riskEXP <- predict(e.CSC, newdata = dt, times = dt$time[25], cause = 1,
                            se = TRUE, iid = TRUE, average.iid = TRUE, product.limit = -1)
    expect_equal(ignore_attr=TRUE,rowMeans(eALL.riskEXP$absRisk.iid[,1,]), eALL.riskEXP$absRisk.average.iid[,1], tol = 1e-6)
    ## check correct  averaging (only R, single obs)
    eALL2.riskEXP <- predict(e.CSC, newdata = dtPred, times = dt$time[25], cause = 1,
                             average.iid = TRUE, product.limit = -1)
    expect_true(max(abs(eALL2.riskEXP$absRisk.average.iid))<1e-6)
    ## check correct  averaging (only R, all)
    eALL2.riskEXP <- predict(e.CSC, newdata = dt, times = dt$time[25], cause = 1,
                             average.iid = TRUE, product.limit = -1)
    expect_equal(ignore_attr=TRUE,eALL.riskEXP$absRisk.average.iid, eALL2.riskEXP$absRisk.average.iid, tol = 1e-6)
    ## check correct  averaging (many timepoints)
    eALL.MriskEXP <- predict(e.CSC, newdata = dt, times = dt$time, cause = 1,
                             se = TRUE, iid = TRUE, average.iid = TRUE, product.limit = -1)
    eALL2.MriskEXP <- predict(e.CSC, newdata = dt, times = dt$time, cause = 1,
                              average.iid = TRUE, product.limit = -1)
    expect_equal(ignore_attr=TRUE,eALL.MriskEXP$absRisk.average.iid, eALL2.MriskEXP$absRisk.average.iid,tol=1e-6)
})
