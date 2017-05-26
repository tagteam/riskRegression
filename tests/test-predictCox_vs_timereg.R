### test-predictCox_vs_timereg.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: maj 26 2017 (17:19) 
##           By: Brice Ozenne
##     Update #: 23
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
    IClambda_GS <- t(as.data.table(m.cox_GS$B.iid))
    m.cox_GSapprox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, rate.sim = 0, clusters = 1:NROW(d), max.timepoint.sim=NULL)
    IClambda_GSapprox <- t(as.data.table(m.cox_GSapprox$B.iid))
  # }}}
  # {{{ no strata, interactions, continous
mI.coxph <- coxph(Surv(eventtime, event) ~ X1*X6, data = d, y = TRUE, x = TRUE)

mI.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(X1) + prop(X6) + prop(X1*X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IClambdaI_GS <- t(as.data.table(mI.cox_GS$B.iid))
  # }}}
  # {{{ no covariate
# cox.aalen do not work without covariable
# m0.cox_GS <- cox.aalen(Surv(eventtime, event) ~ 1, data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
# IClambda0_GS <- t(as.data.table(m0.cox_GS$B.iid))
  
m0.coxph <- coxph(Surv(eventtime, event) ~ 1, data = d, y = TRUE, x = TRUE)
  # }}}
  # {{{ no strata, no interaction, with a categorical variable
mCAT.coxph <- coxph(Surv(eventtime, event) ~ Xcat2 + X6, data = d, y = TRUE, x = TRUE)

mCAT.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(Xcat2) + prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IClambdaCAT_GS <- t(as.data.table(mCAT.cox_GS$B.iid)) 
  # }}}
  # {{{ no strata but with ties 
  mTies.coxph <- coxph(Surv(eventtimeTies, event) ~ X1+X6, data = d3, y = TRUE, x = TRUE)

  mTies.cox_GS <- cox.aalen(Surv(eventtimeTies, event) ~ prop(X1)+prop(X6), data = d3, resample.iid = TRUE, max.timepoint.sim=NULL)
  IClambdaTies_GS <- t(as.data.table(mTies.cox_GS$B.iid))
  # }}}
  # {{{ strata, no interaction, continuous
  mStrata.cox_GS <- cox.aalen(Surv(eventtime, event) ~ strata(St)-1 + prop(X1) + prop(X6), data = dStrata, 
                        resample.iid = TRUE, max.timepoint.sim=NULL)
  
  mStrata.coxph <- coxph(Surv(eventtime, event) ~ strata(St) + X1 + X6, data = dStrata, y = TRUE, x = TRUE)
  # }}}

# }}}


# {{{ 1- Compare iid decomposition (Cox)

  
# {{{ 1a- no strata, no interaction, continous
IC.coxph <- iidCox(m.coxph, keep.times = FALSE)
ICapprox.coxph <- iidCox(m.coxph, keep.times = FALSE, exact = FALSE)
    
IC.coxph_d2 <- iidCox(m.coxph_d2, keep.times = FALSE)
  
IC.cph <- iidCox(m.cph, keep.times = FALSE)
  
test_that("iid beta",{
    expect_equal(unname(IC.coxph$ICbeta),m.cox_GS$gamma.iid)
    expect_equal(unname(IC.coxph$ICbeta),unname(IC.coxph_d2$ICbeta[order(d2$eventtime),]))
    expect_equal(unname(IC.cph$ICbeta),m.cox_GS$gamma.iid, tol = 1e-2)

    expect_equal(unname(ICapprox.coxph$ICbeta),m.cox_GSapprox$gamma.iid)
})
  
test_that("iid lambda0",{
    expect_equal(as.double(IC.coxph$ICcumhazard[[1]]), as.double(IClambda_GS[,-1]))
    expect_equal(as.double(IC.coxph$ICcumhazard[[1]]), as.double(IC.coxph_d2$ICcumhazard[[1]][order(d2$eventtime),]))
    expect_equal(as.double(IC.cph$ICcumhazard[[1]]), as.double(IClambda_GS[,-1]), tol = 1e-4)

    expect_equal(as.double(ICapprox.coxph$ICcumhazard[[1]]), as.double(IClambda_GSapprox[,-1]))
})
# }}}
  
# {{{ 1b- before the first event
data(Melanoma)
Melanoma$status[order(Melanoma$time)]
fitGS <- cox.aalen(Surv(time,status==1)~prop(sex), data=Melanoma)
ICGS_lambda0 <- t(as.data.table(fitGS$B.iid))
  
fit1 <- coxph(Surv(time,status==1)~sex, data=Melanoma, x=TRUE, y=TRUE)
iid1 <- iidCox(fit1)
  
test_that("iid beta - start with censoring",{
    expect_equal(unname(iid1$ICbeta),fitGS$gamma.iid)
})
  
test_that("iid lambda0 - start with censoring",{
    expect_equal(as.double(iid1$ICcumhazard[[1]]), as.double(ICGS_lambda0[,-1]))
})
# }}}

  # {{{ 1c- no strata, interactions, continous  
  ICI.coxph <- iidCox(mI.coxph, keep.times = FALSE)
  
  test_that("iid beta - interaction",{
    expect_equal(unname(ICI.coxph$ICbeta),mI.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - interaction",{
    expect_equal(as.double(ICI.coxph$ICcumhazard[[1]]), as.double(IClambdaI_GS[,-1]))
  })
  
    # }}}

  # {{{ 1d- no covariate
    IC0.coxph <- iidCox(m0.coxph, keep.times = FALSE) # how to test the result
    # }}}
    
  # {{{ 1e- no strata, no interaction, with a categorical variable
  ICCAT.coxph <- iidCox(mCAT.coxph, keep.times = FALSE)
  
  test_that("iid beta - categorical",{
    expect_equal(unname(ICCAT.coxph$ICbeta),mCAT.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - categorical",{
      expect_equal(as.double(ICCAT.coxph$ICcumhazard[[1]]),
                   as.double(IClambdaCAT_GS[,-1]))
  })

    # }}}

  # {{{ 1f- no strata but with ties 
  
  ICTies.coxph <- iidCox(mTies.coxph, keep.times = FALSE)
  
  # test_that("iid beta - categorical",{
  #   expect_equal(unname(ICTies.coxph$ICbeta),mTies.cox_GS$gamma.iid)
  # })
    data.frame(GS = m.cox_GS$gamma.iid, GSTies = mTies.cox_GS$gamma.iid,
               RR = IC.coxph$ICbeta, RRTies = ICTies.coxph$ICbeta)[1:10,]
  
    # test_that("iid lambda0 - categorical",{
    #   expect_equal(as.double(IC.RR$ICcumhazard[[1]]), as.double(IClambda.timereg[,-1]))
    # })
    # }}}

  # {{{ 1g - strata, no interaction, continuous
  
  ICStrata.coxph <- iidCox(mStrata.coxph)
  
  test_that("iid beta - strata",{
    expect_equal(unname(ICStrata.coxph$ICbeta),mStrata.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - strata",{
    
    for(iStrata in 1:length(unique(dStrata$St))){
      IC.GS <- do.call(rbind,
                       lapply(mStrata.cox_GS$B.iid,function(x){x[,iStrata]})
      )
      colnames(IC.GS) <- mStrata.cox_GS$time.sim.resolution
      
      checkTimes <- intersect(mStrata.cox_GS$time.sim.resolution,
                              ICStrata.coxph$time[[iStrata]])
      
      term1 <- ICStrata.coxph$ICcumhazard[[iStrata]][,which(ICStrata.coxph$time[[iStrata]] %in% checkTimes),drop = FALSE]
      term2 <- IC.GS[,which(mStrata.cox_GS$time.sim.resolution %in% checkTimes)]
      diff <- term1 - term2
      expect_true(all(abs(na.omit(as.double(diff)))<1e-10))
    }
    
  })
 
    # }}}


# }}}


# {{{ 3- Compare survival and standard errror of the survival (Cox)
    
# {{{ 3a- no strata, no interaction, continous    
   
test_that("predictionsSE",{
    ## at fix time
    predGS <- predict(m.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(m.coxph, newdata = d, times = 10, se = TRUE)
    predRR2 <- predictCox(m.coxph_d2, newdata = d, times = 10, se = TRUE)

    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))
    
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))

    # at event time
    predGS <- predict(m.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(m.coxph, newdata = d, times = d$eventtime, se = TRUE)
    predRR2 <- predictCox(m.coxph_d2, newdata = d, times = d$eventtime, se = TRUE)

    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR2$survival), as.double(predGS$S0))
    
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    expect_equal(as.double(predRR2$survival.se), as.double(predGS$se.S0))

    # after last event
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e8, se = TRUE)
    expect_true(all(is.na(predRR1$survival)))
    expect_true(all(is.na(predRR1$survival.se)))

    # before first event
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e-8, se = TRUE)
    expect_true(all(predRR1$survival==1))
    expect_true(all(predRR1$survival.se==0))

    # both
    predRR1 <- predictCox(m.coxph, newdata = d, times = c(1e-8,10,1e8), se = TRUE)
    expect_true(all(is.na(predRR1$survival[,3])))
    expect_true(all(predRR1$survival[,1]==1))
    expect_true(all(is.na(predRR1$survival.se[,3])))
    expect_true(all(predRR1$survival.se[,1]==0))
})
  # }}}
  

# {{{ 3b- no strata, interactions, continous  
test_that("predictionsSE - interaction",{
    predGS <- predict(mI.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(mI.coxph, newdata = d, times = 10, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    predGS <- predict(mI.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(mI.coxph, newdata = d, times = d$eventtime, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})
  
# }}}
    
# {{{ 3c- no strata, no interaction, with a categorical variable
test_that("predictionsSE - categorical",{
    predGS <- predict(mCAT.cox_GS, newdata = d, times = 10)
    predRR1 <- predictCox(mCAT.coxph, newdata = d, times = 10, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    predGS <- predict(mCAT.cox_GS, newdata = d, times = d$eventtime)
    predRR1 <- predictCox(mCAT.coxph, newdata = d, times = d$eventtime, se = TRUE)
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})

# }}}

# {{{ 3d- strata, no interaction, continuous

test_that("predictionsSE - strata",{
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 2)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 2, se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
    
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 1:3)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 1:3, se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))

    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = d$eventtime[1:10])
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = d$eventtime[1:10], se = TRUE)
    
    expect_equal(as.double(predRR1$survival), as.double(predGS$S0))
    expect_equal(as.double(predRR1$survival.se), as.double(predGS$se.S0))
})
 
# }}}

# }}}


#----------------------------------------------------------------------
### test-predictCox_vs_timereg.R ends here
