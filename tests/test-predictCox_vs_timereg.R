### test-predictCox_vs_timereg.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: jun  1 2018 (14:42) 
##           By: Brice Ozenne
##     Update #: 75
#----------------------------------------------------------------------
## 
### Commentary: 
## Compare the iid decomposition (lambda and beta) obtained with iidCox and timereg:
## Compare the survival and its standard error obtained with iidCox and timereg:
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(prodlim)
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(timereg)

context("Compare predictCox to timereg")

# {{{ data
set.seed(10)
d <- sampleData(5e1, outcome = "survival")[,.(eventtime,event,X1,X2,X6)]
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
d[ , X16 := X1*X6]
d[ , Xcat2 := as.factor(paste0(X1,X2))]


d2 <- copy(d)

setkey(d,eventtime) # only d is sorted

d3 <- copy(d)[1:10,]
d3[, event := 1]
d3[7:8, X1 := 1]
d3[7:8, X6 := 1]
setkey(d3, eventtime)
d3[, eventtimeTies := eventtime]
d3[7:8, eventtimeTies := eventtime[1]]  # ties

dStrata <- d
dStrata$St <- rbinom(n = NROW(d), size = 2, prob = c(1/3,1/2)) # strata
dStrata2 <- copy(dStrata)
setkeyv(dStrata, c("St", "eventtime"))
# }}}

# {{{ models - Gold standard
# {{{ no strata, no interaction, continous
m.coxph <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
m.coxph_d2 <- coxph(Surv(eventtime, event) ~ X1+X6, data = d2, y = TRUE, x = TRUE)
m.cph <- cph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
    
m.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IFlambda_GS <- t(as.data.table(m.cox_GS$B.iid))
m.cox_GSapprox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, rate.sim = 0, clusters = 1:NROW(d), max.timepoint.sim=NULL)
IFlambda_GSapprox <- t(as.data.table(m.cox_GSapprox$B.iid))
# }}}
# {{{ no strata, interactions, continous
mI.coxph <- coxph(Surv(eventtime, event) ~ X1*X6, data = d, y = TRUE, x = TRUE)
mI.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(X1) + prop(X6) + prop(X1*X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IFlambdaI_GS <- t(as.data.table(mI.cox_GS$B.iid))
# }}}

# {{{ no covariate
# cox.aalen do not work without covariable
# m0.cox_GS <- cox.aalen(Surv(eventtime, event) ~ 1, data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
# IFlambda0_GS <- t(as.data.table(m0.cox_GS$B.iid))
m0.coxph <- coxph(Surv(eventtime, event) ~ 1, data = d, y = TRUE, x = TRUE)
# }}}
# {{{ no strata, no interaction, with a categorical variable
mCAT.coxph <- coxph(Surv(eventtime, event) ~ Xcat2 + X6, data = d, y = TRUE, x = TRUE)

mCAT.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(Xcat2) + prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IFlambdaCAT_GS <- t(as.data.table(mCAT.cox_GS$B.iid)) 
  # }}}
  # {{{ no strata but with ties 
  mTies.coxph <- coxph(Surv(eventtimeTies, event) ~ X1+X6, data = d3, y = TRUE, x = TRUE)

  mTies.cox_GS <- cox.aalen(Surv(eventtimeTies, event) ~ prop(X1)+prop(X6), data = d3, resample.iid = TRUE, max.timepoint.sim=NULL)
  IFlambdaTies_GS <- t(as.data.table(mTies.cox_GS$B.iid))
  # }}}
  # {{{ strata, no interaction, continuous
  mStrata.cox_GS <- cox.aalen(Surv(eventtime, event) ~ strata(St)-1 + prop(X1) + prop(X6), data = dStrata, 
                        resample.iid = TRUE, max.timepoint.sim=NULL)
  
  mStrata.coxph <- coxph(Surv(eventtime, event) ~ strata(St) + X1 + X6, data = dStrata, y = TRUE, x = TRUE)
  # }}}

# }}}

# {{{ 1- Compare iid decomposition (Cox)

  
# {{{ 1a- no strata, no interaction, continous
IF.coxph <- iidCox(m.coxph, keep.times = FALSE)
IFapprox.coxph <- iidCox(m.coxph, keep.times = FALSE, store.iid = "approx")
IFminimal.coxph <- iidCox(m.coxph, keep.times = FALSE, store.iid = "minimal")
IF.coxph_d2 <- iidCox(m.coxph_d2, keep.times = FALSE)
IF.cph <- iidCox(m.cph, keep.times = FALSE)
test_that("iid beta",{
    expect_equal(unname(IF.coxph$IFbeta),m.cox_GS$gamma.iid)
    expect_equal(unname(IF.coxph$IFbeta),unname(IF.coxph_d2$IFbeta[order(d2$eventtime),]))
    expect_equal(unname(IF.cph$IFbeta),m.cox_GS$gamma.iid, tol = 1e-2)
    expect_equal(unname(IFapprox.coxph$IFbeta),m.cox_GSapprox$gamma.iid)
    expect_equal(unname(IFminimal.coxph$IFbeta),m.cox_GS$gamma.iid)
})
  
test_that("iid lambda0",{
    expect_equal(as.double(IF.coxph$IFcumhazard[[1]]), as.double(IFlambda_GS[,-1]))
    expect_equal(as.double(IF.coxph$IFcumhazard[[1]]), as.double(IF.coxph_d2$IFcumhazard[[1]][order(d2$eventtime),]))
    expect_equal(as.double(IF.cph$IFcumhazard[[1]]), as.double(IFlambda_GS[,-1]), tol = 1e-4)
    expect_equal(as.double(IFapprox.coxph$IFcumhazard[[1]]), as.double(IFlambda_GSapprox[,-1]))
})
# }}}
  
# {{{ 1b- before the first event
data(Melanoma)
# Melanoma$status[order(Melanoma$time)]
fitGS <- cox.aalen(Surv(time,status==1)~prop(sex), data=Melanoma)
IFGS_lambda0 <- t(as.data.table(fitGS$B.iid))
fit1 <- coxph(Surv(time,status==1)~sex, data=Melanoma, x=TRUE, y=TRUE)
iid1 <- iidCox(fit1)
test_that("iid beta - start with censoring",{
    expect_equal(unname(iid1$IFbeta),fitGS$gamma.iid)
})
test_that("iid lambda0 - start with censoring",{
    expect_equal(as.double(iid1$IFcumhazard[[1]]), as.double(IFGS_lambda0[,-1]))
})
# }}}

  # {{{ 1c- no strata, interactions, continous  
  IFI.coxph <- iidCox(mI.coxph, keep.times = FALSE)
  test_that("iid beta - interaction",{
    expect_equal(unname(IFI.coxph$IFbeta),mI.cox_GS$gamma.iid)
  })
  test_that("iid lambda0 - interaction",{
    expect_equal(as.double(IFI.coxph$IFcumhazard[[1]]), as.double(IFlambdaI_GS[,-1]))
  })
    # }}}

# {{{ 1d- no covariate
IF0.coxph <- iidCox(m0.coxph, keep.times = FALSE) # how to test the result
# }}}
    
# {{{ 1e- no strata, no interaction, with a categorical variable
IFCAT.coxph <- iidCox(mCAT.coxph, keep.times = FALSE)
test_that("iid beta - categorical",{
    expect_equal(unname(IFCAT.coxph$IFbeta),mCAT.cox_GS$gamma.iid)
})
  
test_that("iid lambda0 - categorical",{
    expect_equal(as.double(IFCAT.coxph$IFcumhazard[[1]]),
                 as.double(IFlambdaCAT_GS[,-1]))
})

    # }}}

  # {{{ 1f- no strata but with ties 
  
  IFTies.coxph <- iidCox(mTies.coxph, keep.times = FALSE)
  
# test_that("iid beta - categorical",{
#   expect_equal(unname(IFTies.coxph$IFbeta),mTies.cox_GS$gamma.iid)
# })
## data.frame(GS = m.cox_GS$gamma.iid, GSTies = mTies.cox_GS$gamma.iid,
## RR = IF.coxph$IFbeta, RRTies = IFTies.coxph$IFbeta)[1:10,]
  
    # test_that("iid lambda0 - categorical",{
    #   expect_equal(as.double(IF.RR$IFcumhazard[[1]]), as.double(IFlambda.timereg[,-1]))
    # })
    # }}}

  # {{{ 1g - strata, no interaction, continuous
  
  IFStrata.coxph <- iidCox(mStrata.coxph)
  
  test_that("iid beta - strata",{
    expect_equal(unname(IFStrata.coxph$IFbeta),mStrata.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - strata",{
    
    for(iStrata in 1:length(unique(dStrata$St))){
      IF.GS <- do.call(rbind,
                       lapply(mStrata.cox_GS$B.iid,function(x){x[,iStrata]})
      )
      colnames(IF.GS) <- mStrata.cox_GS$time.sim.resolution
      
      checkTimes <- intersect(mStrata.cox_GS$time.sim.resolution,
                              IFStrata.coxph$time[[iStrata]])
      
      term1 <- IFStrata.coxph$IFcumhazard[[iStrata]][,which(IFStrata.coxph$time[[iStrata]] %in% checkTimes),drop = FALSE]
      term2 <- IF.GS[,which(mStrata.cox_GS$time.sim.resolution %in% checkTimes)]
      diff <- term1 - term2
      expect_true(all(abs(na.omit(as.double(diff)))<1e-10))
    }
    
  })
 
    # }}}


# }}}

# {{{ 2- Compare survival and standard errror of the survival (Cox)
    
# {{{ 2a- no strata, no interaction, continous    
   
test_that("predictionsSE",{

#### at fix time
    predGS <- predict(m.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(m.coxph, newdata = d, times = 10, se = TRUE, band = TRUE)
    predRR2 <- predictCox(m.coxph_d2, newdata = d, times = 10, se = TRUE, band = TRUE)

    predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10)
    predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10)
    
    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.101588, 0.070717), 
                     "survival.se" = c(0.086529, 0.057379), 
                     "survival.lower" = c(0, 0), 
                     "survival.upper" = c(0.271182, 0.183179),                     
                     "survival.quantileBand" = c(1.963483, 1.982114),
                     "survival.lowerBand" = c(0, 0), 
                     "survival.upperBand" = c(0.271486, 0.18445))
    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.10159, 0.07072), 
                     "survival.se" = c(0.37246, 0.30629), 
                     "survival.lower" = c(0.00869, 0.008), 
                     "survival.upper" = c(0.3322, 0.23378),
                     "survival.quantileBand" = c(1.96348, 1.98211),
                     "survival.lowerBand" = c(0.00864, 0.00774), 
                     "survival.upperBand" = c(0.33268, 0.23609))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE],
                 GS, tol = 1e-4, scale = 1)

    
#### at event time
    predGS <- predict(m.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(m.coxph, newdata = d, times = d$eventtime, se = TRUE, band = TRUE)
    predRR2 <- predictCox(m.coxph_d2, newdata = d, times = d$eventtime, se = TRUE, band = TRUE)

    predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10)
    predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10)

    ## punctual estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))

    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.170031, 0.170031), 
                     "survival" = c(0.983427, 0.980827), 
                     "survival.se" = c(0.016414, 0.020011), 
                     "survival.lower" = c(0.951257, 0.941606), 
                     "survival.upper" = c(1, 1), 
                     "survival.quantileBand" = c(2.70855, 2.896028), 
                     "survival.lowerBand" = c(0.93897, 0.922875), 
                     "survival.upperBand" = c(1, 1))
    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.17003, 0.17003), 
                     "survival" = c(0.98343, 0.98083), 
                     "survival.se" = c(0.99872, 1.0539), 
                     "survival.lower" = c(0.88839, 0.85835), 
                     "survival.upper" = c(0.99764, 0.99755), 
                     "survival.quantileBand" = c(2.70855, 2.89603), 
                     "survival.lowerBand" = c(0.77885, 0.66389), 
                     "survival.upperBand" = c(0.99888, 0.99909))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

#### after last event
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e8, se = TRUE, band = TRUE)
    expect_true(all(is.na(predRR1$survival)))
    expect_true(all(is.na(predRR1$survival.se)))
    
#### before first event
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e-8, se = TRUE)
    expect_true(all(predRR1$survival==1))
    expect_true(all(predRR1$survival.se==0))

#### both
    predRR1 <- predictCox(m.coxph, newdata = d, times = c(1e-8,10,1e8), se = TRUE)
    expect_true(all(is.na(predRR1$survival[,3])))
    expect_true(all(predRR1$survival[,1]==1))
    expect_true(all(is.na(predRR1$survival.se[,3])))
    expect_true(all(predRR1$survival.se[,1]==0))
})
                                        # }}}
  

# {{{ 2b- no strata, interactions, continous  
test_that("predictionsSE - interaction",{

#### at fix time
    predGS <- predict(mI.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(mI.coxph, newdata = d, times = 10, se = TRUE, band = TRUE)

    predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10)
    predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10)

    ## punctial estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    
    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.113587, 0.075708), 
                     "survival.se" = c(0.102236, 0.061297), 
                     "survival.lower" = c(0, 0), 
                     "survival.upper" = c(0.313966, 0.195847), 
                     "survival.quantileBand" = c(1.950346, 1.960436), 
                     "survival.lowerBand" = c(0, 0), 
                     "survival.upperBand" = c(0.312983, 0.195876))
    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(10, 10), 
                     "survival" = c(0.11359, 0.07571), 
                     "survival.se" = c(0.41379, 0.31371), 
                     "survival.lower" = c(0.00749, 0.00845), 
                     "survival.upper" = c(0.38035, 0.2477), 
                     "survival.quantileBand" = c(1.95035, 1.96044), 
                     "survival.lowerBand" = c(0.00763, 0.00845), 
                     "survival.upperBand" = c(0.37888, 0.24776))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

#### at event time
    predGS <- predict(mI.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(mI.coxph, newdata = d, times = d$eventtime, se = TRUE, band = TRUE)

    predRR1.none <- confint(predRR1, survival.transform = "none", seed = 10)
    predRR1.loglog <- confint(predRR1, survival.transform = "loglog", seed = 10)

    ## punctial estimate
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    
    ## standard error
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    ## confidence intervals/band
    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.170031, 0.170031), 
                     "survival" = c(0.985361, 0.982654), 
                     "survival.se" = c(0.016408, 0.019694), 
                     "survival.lower" = c(0.953201, 0.944054), 
                     "survival.upper" = c(1, 1), 
                     "survival.quantileBand" = c(2.655669, 2.846195), 
                     "survival.lowerBand" = c(0.941785, 0.926601), 
                     "survival.upperBand" = c(1, 1))
    ## butils::object2script(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], digit = 6)
    expect_equal(as.data.table(predRR1.none)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)

    GS <- data.table("observation" = c(6, 7), 
                     "times" = c(0.17003, 0.17003), 
                     "survival" = c(0.98536, 0.98265), 
                     "survival.se" = c(1.12917, 1.14539), 
                     "survival.lower" = c(0.87384, 0.84775), 
                     "survival.upper" = c(0.99839, 0.99815), 
                     "survival.quantileBand" = c(2.65567, 2.84619), 
                     "survival.lowerBand" = c(0.74392, 0.63393), 
                     "survival.upperBand" = c(0.99927, 0.99933))
    ## butils::object2script(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], digit = 5)
    expect_equal(as.data.table(predRR1.loglog)[6:7,names(GS),with=FALSE], GS, tol = 1e-4, scale = 1)
})
  
# }}}
    
# {{{ 2c- no strata, no interaction, with a categorical variable
test_that("predictionsSE - categorical",{

#### at fix time
    predGS <- predict(mCAT.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(mCAT.coxph, newdata = d, times = 10, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

#### at event time
    predGS <- predict(mCAT.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(mCAT.coxph, newdata = d, times = d$eventtime, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})

# }}}

# {{{ 3d- strata, no interaction, continuous

test_that("predictionsSE - strata",{

#### at fix time
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 2)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 2, se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 1:3)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 1:3, se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

#### at event time
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = d$eventtime[1:10])
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = d$eventtime[1:10], se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})
 
# }}}

# }}}

# {{{ 3- IChazard vs ICcumhazard

IF.coxph <- iidCox(m.coxph, keep.times = FALSE)

expect_equal(IF.coxph$IFcumhazard[[1]],t(apply(IF.coxph$IFhazard[[1]],1,cumsum)))


# }}}

# {{{ 4- Fast computation standard error
set.seed(10)
d <- sampleData(50, outcome = "survival")
setkey(d,time)

m.coxph <- coxph(Surv(time, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

## system.time(
##     res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata, store.iid = "minimal", se = TRUE, iid = FALSE)
## )
## system.time(
##     res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata, store.iid = "full", se = TRUE, iid = FALSE)
## )

test_that("iid minimal - no strata", {
    res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", se = TRUE, iid = TRUE) 
    res2 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", average.iid = TRUE) 
    res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$cumhazard.se,res3$cumhazard.se)
    expect_equal(res1$survival.se,res3$survival.se)
    expect_equal(res1$cumhazard.iid,res3$cumhazard.iid)
    expect_equal(res1$survival.iid,res3$survival.iid)

    expect_equal(res2$cumhazard.average.iid, t(apply(res3$cumhazard.iid,2:3,mean)))
    expect_equal(res2$survival.average.iid, t(apply(res3$survival.iid,2:3,mean)))
})

m.coxph <- coxph(Surv(time, event) ~ strata(X1)+X6, data = d, y = TRUE, x = TRUE)
 
seqTime <- c(1e-16,4:10,d$time[1:10],1e6)
newdata <- d

test_that("iid minimal - strata", {
    newdata <- rbind(d[1],d[1])
    res1 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", se = TRUE, iid = TRUE) 
    res2 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "minimal", average.iid = TRUE) 
    res3 <- predictCox(m.coxph, times = seqTime, newdata = newdata,
                       type = c("cumhazard", "survival"),
                       store.iid = "full", se = TRUE, iid = TRUE)
    expect_equal(res1$cumhazard.se,res3$cumhazard.se)
    expect_equal(res1$survival.se,res3$survival.se)
    expect_equal(res1$cumhazard.iid,res3$cumhazard.iid)
    expect_equal(res1$survival.iid,res3$survival.iid)

    expect_equal(res2$cumhazard.average.iid, t(apply(res3$cumhazard.iid,2:3,mean)))
    expect_equal(res2$survival.average.iid, t(apply(res3$survival.iid,2:3,mean)))
})

# }}}

# {{{ 5- Time varying variables
# NOTE: I am not sure what timereg is exactly doing in this case

## d <- sampleData(1e2, outcome = "survival")
## d$start <- runif(NROW(d),min=0,max=(d$eventtime-0.1) )

## test_that("predicted survival with time varying variables",{
    ## fit.coxph <- coxph(Surv(start, eventtime, event) ~ X1, data = d, x = TRUE, y = TRUE)
    ## fit.timereg <- cox.aalen(Surv(start, eventtime, event) ~ prop(X1), data = d)
    ## predTimes <- sort(unique(d$eventtime))
    ## M1 <- predictCox(fit.coxph, newdata = d, time = predTimes)$survival
    ## M2 <- predict(fit.timereg, newdata = d, time = predTimes)$S0
    ## expect_equal(M1,M2)
## })
# }}}
#----------------------------------------------------------------------
### test-predictCox_vs_timereg.R ends here

