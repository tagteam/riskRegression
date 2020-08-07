### test-iidCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: aug  7 2020 (12:24) 
##           By: Brice Ozenne
##     Update #: 117
#----------------------------------------------------------------------
## 
### Commentary: 
## Check internal consistency of the iid decomposition (e.g. after last event time, IF cumhazard vs. IF hazard)
## Compare the iid decomposition (lambda and beta) obtained with iidCox and timereg
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Settings
library(testthat)
library(riskRegression)
library(survival)
library(rms)
library(timereg)
library(data.table)

context("function iidCox")

## * Internal consistency
cat("[iidCox] internal consistency \n")

## ** Data
set.seed(10)
dt <- sampleData(20, outcome = "survival")[,.(time,event,X1,X2)]

timeFirstEvent <- dt[event==1,min(time)]
timeLastEvent <- dt[event==1,max(time)]
timeLastObs <- dt[,max(time)]
vec.times <-  timeFirstEvent+1:5
dt.newdata <- dt[1:4]

## ** Model
coxph.fit <- coxph(Surv(time,event==1)~ X1+X2,
                   data=dt, x = TRUE, y = TRUE)
cph.fit <- coxph(Surv(time,event==1)~ X1+X2,
                 data=dt, x = TRUE, y = TRUE)

## ** Tests
iid.coxph <-  iidCox(coxph.fit, return.object = FALSE)

test_that("[iidCox] consistency coxph,cph",{
    iid.cph <-  iidCox(cph.fit, return.object = FALSE)
    expect_equal(iid.cph, iid.coxph)
})

test_that("[iidCox] consistency with manual specification of newdata and times",{
    iid.test <- iidCox(coxph.fit, newdata = dt, tau.hazard = iid.coxph$time[[1]], return.object = FALSE)
    expect_equal(iid.coxph, iid.test)

    new.order <- sample.int(length(iid.coxph$time[[1]]))
    iid.test <- iidCox(coxph.fit, newdata = dt, tau.hazard = iid.coxph$time[[1]][new.order], return.object = FALSE)
    expect_equal(iid.coxph$IFhazard[[1]], iid.test$IFhazard[[1]][,colnames(iid.coxph$IFhazard[[1]])])
    expect_equal(iid.coxph$IFcumhazard[[1]], iid.test$IFcumhazard[[1]][,colnames(iid.coxph$IFcumhazard[[1]])])
    expect_equal(iid.coxph$IFbeta, iid.coxph$IFbeta)
})

test_that("[iidCox] cumsum(iid hazard) = iid cumhazard",{
    M.test <- t(apply(iid.coxph$IFhazard[[1]], 1, cumsum))
    expect_equal(M.test, iid.coxph$IFcumhazard[[1]])
})

test_that("[iidCox] iid hazard = 0 at non event times and NA after last observation",{
    iid.test <- iidCox(coxph.fit, tau.hazard = sort(c(dt$time,c(0.1,1,5,8,12,25))), return.object = FALSE)

    expect_equal(iid.test$IFhazard[[1]][,as.character(iid.coxph$time[[1]])],
                 iid.coxph$IFhazard[[1]])

    vec.noEventTime <- setdiff(iid.test$time[[1]],dt[event==1,time])
    vec.noEventTime_beforeTmax <- vec.noEventTime[vec.noEventTime<=timeLastObs]
    vec.noEventTime_afterTmax <- vec.noEventTime[vec.noEventTime>timeLastObs]
    
    expect_true(all(iid.test$IFhazard[[1]][,as.character(vec.noEventTime_beforeTmax)]==0))
    expect_true(all(is.na(iid.test$IFhazard[[1]][,as.character(vec.noEventTime_afterTmax)])))

    ## additional tests
    ## compatibility IFhazard, IFcumhazard
    M.test <- t(apply(iid.test$IFhazard[[1]], 1, cumsum))
    expect_equal(M.test, iid.test$IFcumhazard[[1]])

    ## iid cumHazard[lastEvent+] =  cumHazard[lastEvent]
    vec.time <- intersect(iid.test$time[[1]][iid.test$time[[1]]>timeLastEvent],
                          iid.test$time[[1]][iid.test$time[[1]]<=timeLastObs])
    expect_true(all(iid.test$IFcumhazard[[1]][,as.character(vec.time)]-iid.test$IFcumhazard[[1]][,as.character(timeLastEvent)]==0))
})
## ** selectJump
## test_that("[iidCox] selectJump",{
##     seqTest <- c(0,dt$time[1:10],dt$time[1:10]-1e-12,dt$time[1:10]+1e-5,1e5,1)
##     GS <-  iidCox(coxph.fit,
##                   tau.hazard = seqTest,
##                   return.object = FALSE)
##     test <- riskRegression:::selectJump(iid.coxph, times = seqTest, type = c("hazard","cumhazard"))

##     expect_equal(unname(GS$IFhazard[[1]]),unname(test$IFhazard[[1]]))
##     expect_equal(unname(GS$IFcumhazard[[1]]),unname(test$IFcumhazard[[1]]))
## })

## ** Empty strata
set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

e.coxph <- coxph(Surv(time, event) ~ strata(X1), data = dtS,
                 x = TRUE, y = TRUE)
test_that("[iidCox] empty strata", {
    test <- iidCox(e.coxph, tau.max = 1)
    expect_equal(test$iid$time[[1]],0)
    expect_equal(test$iid$time[[2]],c(0.2,0.6,0.9))
    expect_equal(test$iid$time[[3]],0)

    expect_equal(unname(test$iid$IFhazard[[1]]), matrix(0, nrow = NROW(dtS), ncol = 1))
    GS <- cbind("0.2" = c(0, 0, 0, -0.00173611, 0.03993056, 0, 0, -0.00173611, 0, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611, 0, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611, 0, -0.00173611, -0.00173611, -0.00173611, -0.00173611, 0, 0, -0.00173611, 0, 0, -0.00173611, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611), 
                "0.6" = c(0, 0, 0, -0.00189036, 0, 0, 0, -0.00189036, 0, 0, 0, 0, -0.00189036, -0.00189036, -0.00189036, 0, -0.00189036, 0, 0, 0, 0, -0.00189036, -0.00189036, -0.00189036, 0, 0, 0, -0.00189036, -0.00189036, -0.00189036, 0, -0.00189036, 0, -0.00189036, -0.00189036, -0.00189036, -0.00189036, 0, 0, -0.00189036, 0, 0, -0.00189036, 0, 0, -0.00189036, 0.0415879, -0.00189036, 0, -0.00189036), 
                "0.9" = c(0, 0, 0, -0.00206612, 0, 0, 0, -0.00206612, 0, 0, 0, 0, -0.00206612, -0.00206612, -0.00206612, 0, -0.00206612, 0, 0, 0, 0, 0.04338843, -0.00206612, -0.00206612, 0, 0, 0, -0.00206612, -0.00206612, -0.00206612, 0, -0.00206612, 0, -0.00206612, -0.00206612, -0.00206612, -0.00206612, 0, 0, -0.00206612, 0, 0, -0.00206612, 0, 0, -0.00206612, 0, -0.00206612, 0, -0.00206612))
    expect_equal(test$iid$IFhazard[[2]], GS, tol = 1e-6)
    expect_equal(unname(test$iid$IFhazard[[3]]), matrix(0, nrow = NROW(dtS), ncol = 1))

    expect_equal(unname(test$iid$IFcumhazard[[1]]), matrix(0, nrow = NROW(dtS), ncol = 1))
    GS <- cbind("0.2" = c(0, 0, 0, -0.00173611, 0.03993056, 0, 0, -0.00173611, 0, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611, 0, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611, 0, -0.00173611, -0.00173611, -0.00173611, -0.00173611, 0, 0, -0.00173611, 0, 0, -0.00173611, 0, 0, -0.00173611, -0.00173611, -0.00173611, 0, -0.00173611), 
                "0.6" = c(0, 0, 0, -0.00362647, 0.03993056, 0, 0, -0.00362647, 0, 0, 0, 0, -0.00362647, -0.00362647, -0.00362647, 0, -0.00362647, 0, 0, 0, 0, -0.00362647, -0.00362647, -0.00362647, 0, 0, 0, -0.00362647, -0.00362647, -0.00362647, 0, -0.00362647, 0, -0.00362647, -0.00362647, -0.00362647, -0.00362647, 0, 0, -0.00362647, 0, 0, -0.00362647, 0, 0, -0.00362647, 0.03985179, -0.00362647, 0, -0.00362647), 
                "0.9" = c(0, 0, 0, -0.00569259, 0.03993056, 0, 0, -0.00569259, 0, 0, 0, 0, -0.00569259, -0.00569259, -0.00569259, 0, -0.00569259, 0, 0, 0, 0, 0.03976196, -0.00569259, -0.00569259, 0, 0, 0, -0.00569259, -0.00569259, -0.00569259, 0, -0.00569259, 0, -0.00569259, -0.00569259, -0.00569259, -0.00569259, 0, 0, -0.00569259, 0, 0, -0.00569259, 0, 0, -0.00569259, 0.03985179, -0.00569259, 0, -0.00569259))
    expect_equal(test$iid$IFcumhazard[[2]], GS, tol = 1e-6)
    expect_equal(unname(test$iid$IFcumhazard[[3]]), matrix(0, nrow = NROW(dtS), ncol = 1))

    expect_equal(test$iid$etime.max, c(11.4,14.8,9.2))
    ## expect_equal(test$iid$indexObs,order(dtS$time))
})

## * Compare to timereg
## ** Data
data(Melanoma)

set.seed(10)
dt <- sampleData(5e1, outcome = "survival")[,.(time,event,X1,X2,X6)]
dt[,X1:=as.numeric(as.character(X1))]
dt[,X2:=as.numeric(as.character(X2))]
dt[ , X16 := X1*X6]
dt[ , Xcat2 := as.factor(paste0(X1,X2))]

## sorted dataset
dt.sort <- copy(dt)
setkeyv(dt.sort,c("time")) 

## dataset with ties
dtTies <- copy(dt)[1:10,]
dtTies[, event := 1]
dtTies[7:8, X1 := 1]
dtTies[7:8, X6 := 1]
setkeyv(dtTies, c("time"))
dtTies[, timeTies := time]
dtTies[7:8, timeTies := time[1]]  # ties

dtStrata <- copy(dt)
dtStrata[, strata :=  rbinom(n = .N, size = 2, prob = c(1/3,1/2))] # strata
dtStrata.sort <- copy(dtStrata)
setkeyv(dtStrata.sort, c("strata", "time"))

## ** No strata, no interaction, continous
cat("[iidCox] compare to timereg - no strata, no interaction, continuous \n")
## *** Models
e.coxph <- coxph(Surv(time, event) ~ X1+X6, data = dt, y = TRUE, x = TRUE)
e.coxph_sort <- coxph(Surv(time, event) ~ X1+X6, data = dt.sort, y = TRUE, x = TRUE)
    
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1)+prop(X6),
                          data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

eApprox.timereg <- cox.aalen(Surv(time, event) ~ prop(X1)+prop(X6),
                                data = dt, resample.iid = TRUE, max.timepoint.sim=NULL,
                                rate.sim = 0, clusters = 1:NROW(dt))

## *** Extract information
## compute the IF of the baseline hazard and then for the survival
IF.coxph <- iidCox(e.coxph, keep.times = FALSE, store.iid = "full", return.object = FALSE)
IF.coxph_sort <- iidCox(e.coxph_sort, keep.times = FALSE, store.iid = "full", return.object = FALSE)
## store "sufficient statistics" for the IF of the baseline hazard
## and then compute the IF for the survival
IFminimal.coxph <- iidCox(e.coxph, keep.times = FALSE, store.iid = "minimal", return.object = FALSE)
## same as before but using an approximation

IFlambda_GS <- t(as.data.table(e.timereg$B.iid))
 
## *** Tests
test_that("[iidCox] beta - no strata, no interaction, continuous",{
    expect_equal(unname(IF.coxph$IFbeta),e.timereg$gamma.iid)
    expect_equal(unname(IF.coxph$IFbeta[order(dt$time),]),unname(IF.coxph_sort$IFbeta))
    expect_equal(unname(IFminimal.coxph$IFbeta),e.timereg$gamma.iid)
})

test_that("[iidCox] lambda - no strata, no interaction, continuous",{
    expect_equal(as.double(IF.coxph$IFcumhazard[[1]]), as.double(IFlambda_GS[,-1]))
    expect_equal(IF.coxph$IFcumhazard[[1]][order(dt$time),], IF.coxph_sort$IFcumhazard[[1]])
    expect_equal(IFminimal.coxph$IFcumhazard[[1]], NULL)
})
  

## ** no strata, interactions, continous
cat("[iidCox] compare to timereg - no strata, interaction, continuous \n")

## *** Models
e.coxph <- coxph(Surv(time, event) ~ X1*X6, data = dt, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X6) + prop(X1*X6),
                       data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

## *** Extract information
IF.coxph <- iidCox(e.coxph, keep.times = FALSE, return.object = FALSE)
IFlambda_GS <- t(as.data.table(e.timereg$B.iid))

## *** Tests
test_that("[iidCox] beta - no strata, interactions, continuous",{
    expect_equal(unname(IF.coxph$IFbeta),e.timereg$gamma.iid)
})
test_that("[iidCox] lambda - no strata, interactions, continuous",{
    expect_equal(as.double(IF.coxph$IFcumhazard[[1]]), as.double(IFlambda_GS[,-1]))
})


## ** no covariate
cat("[iidCox] compare to timereg - no covariate \n")
## cox.aalen do not work without covariable
## cox.aalen(Surv(eventtime, event) ~ 1, data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)
## Error in model.frame.default(formula = Surv(survs$stop, survs$status) ~  : 
##   type (NULL) incorrect pour la variable 'Z'

test_that("[iidCox] beta - no covariate",{
    e.coxph <- coxph(Surv(time, event) ~ 1, data = dt, y = TRUE, x = TRUE)
    IF.coxph <- iidCox(e.coxph, keep.times = FALSE, return.object = FALSE)
    expect_true(all(is.na(IF.coxph$IFbeta)))
})

## ** no strata, no interaction, with a categorical variable
cat("[iidCox] compare to timereg - no strata, no interaction, categorical \n")

## *** Models
e.coxph <- coxph(Surv(time, event) ~ Xcat2 + X6, data = dt, y = TRUE, x = TRUE)
e.timereg <- cox.aalen(Surv(time, event) ~ prop(Xcat2) + prop(X6), data = dt, resample.iid = TRUE, max.timepoint.sim=NULL)

## coef(e.coxph)
## coef(e.timereg)

## *** Extract information
IFlambda_GS <- t(as.data.table(e.timereg$B.iid))
IF.coxph <- iidCox(e.coxph, keep.times = FALSE, return.object = FALSE)

## *** Tests
test_that("[iidCox] beta - no strata, interactions, categorical",{
    expect_equal(unname(IF.coxph$IFbeta),e.timereg$gamma.iid)
})
  
test_that("[iidCox] lambda - no strata, interactions, categorical",{
    expect_equal(as.double(IF.coxph$IFcumhazard[[1]]),
                 as.double(IFlambda_GS[,-1]))
})



## ** Ties
cat("[iidCox] compare to timereg - ties \n")

e.coxph_Efron <- coxph(Surv(timeTies, event) ~ X1+X6, data = dtTies, y = TRUE, x = TRUE,
                 ties = "efron")
e.coxph_Breslow <- coxph(Surv(timeTies, event) ~ X1+X6, data = dtTies, y = TRUE, x = TRUE,
                         ties = "breslow")
e.timereg <- cox.aalen(Surv(timeTies, event) ~ prop(X1)+prop(X6), data = dtTies, resample.iid = TRUE, max.timepoint.sim=NULL)

## different coef for timereg - cannot compare
coef(e.coxph_Breslow)
coef(e.coxph_Efron)
coef(e.timereg)[,"Coef."]

## *** Extract information
IF.coxph_Efron <- iidCox(e.coxph_Efron, keep.times = FALSE, return.object = FALSE)
IF.coxph_Breslow <- iidCox(e.coxph_Breslow, keep.times = FALSE, return.object = FALSE)

## *** Tests
## ??? Gold standard


## ** Strata, no interaction, continuous
cat("[iidCox] compare to timereg - strata \n")
## *** Model
e.timereg <- cox.aalen(Surv(time, event) ~ strata(strata)-1 + prop(X1) + prop(X6), data = dtStrata, 
                       resample.iid = TRUE, max.timepoint.sim=NULL)
e.coxph <- coxph(Surv(time, event) ~ strata(strata) + X1 + X6, data = dtStrata, y = TRUE, x = TRUE)

## *** Extract information
IF.coxph <- iidCox(e.coxph, return.object = FALSE)

## *** Tests
test_that("[iidCox] beta - strata",{
    expect_equal(unname(IF.coxph$IFbeta),e.timereg$gamma.iid)
})
  
test_that("[iidCox] lambda - strata",{

    name.strata <- unique(dtStrata$strata)
    n.strata <- length(name.strata)
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        IF.GS <- do.call(rbind,
                         lapply(e.timereg$B.iid,function(x){x[,iStrata]})
                         )
        colnames(IF.GS) <- e.timereg$time.sim.resolution
      
        checkTimes <- intersect(e.timereg$time.sim.resolution,
                                IF.coxph$time[[iStrata]])
      
        term1 <- IF.coxph$IFcumhazard[[iStrata]][,which(IF.coxph$time[[iStrata]] %in% checkTimes),drop = FALSE]
        term2 <- IF.GS[,which(e.timereg$time.sim.resolution %in% checkTimes)]
        diff <- term1 - term2
        expect_true(all(abs(na.omit(as.double(diff)))<1e-10))
    }
    
  })



## ** Melanoma data
cat("[iidCox] compare to timereg - Melanoma \n")
data(Melanoma, package = "riskRegression")

test_that("[iidCox] Compare to timereg on Melanoma dta",{
    e.timereg <- cox.aalen(Surv(time,status==1)~prop(sex), data = Melanoma)
    e.coxph <- coxph(Surv(time,status==1)~sex, data=Melanoma, x=TRUE, y=TRUE)

    RR.iid <- iidCox(e.coxph, return.object = FALSE)
    timereg.iidLambda <- t(as.data.table(e.timereg$B.iid))

    expect_equal(unname(RR.iid$IFbeta),e.timereg$gamma.iid)
    expect_equal(as.double(RR.iid$IFcumhazard[[1]]),
                 as.double(timereg.iidLambda[,-1]))
})

#----------------------------------------------------------------------
### test-iidCox.R ends here
