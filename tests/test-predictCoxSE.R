library(prodlim)
library(riskRegression)
library(testthat)
library(rms)
library(survival)

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
if(require(timereg)){
  # {{{ no strata, no interaction, continous
m.coxph <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
m.coxph_d2 <- coxph(Surv(eventtime, event) ~ X1+X6, data = d2, y = TRUE, x = TRUE)
m.cph <- cph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
  
m.cox_GS <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
IClambda_GS <- t(as.data.table(m.cox_GS$B.iid))
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
}
# }}}


# {{{ Cox influence function

if(require(timereg)){
  
    # {{{ no strata, no interaction, continous
  IC.coxph <- iidCox(m.coxph, keep.times = FALSE)
    
  IC.coxph_d2 <- iidCox(m.coxph_d2, keep.times = FALSE)
  
  IC.cph <- iidCox(m.cph, keep.times = FALSE)
  
  test_that("iid beta",{
    expect_equal(unname(IC.coxph$ICbeta),m.cox_GS$gamma.iid)
    expect_equal(unname(IC.coxph$ICbeta),unname(IC.coxph_d2$ICbeta[order(d2$eventtime),]))
    expect_equal(unname(IC.cph$ICbeta),m.cox_GS$gamma.iid, tol = 1e-2)
  })
  
    test_that("iid lambda0",{
    expect_equal(as.double(IC.coxph$ICcumhazard[[1]]), as.double(IClambda_GS[,-1]))
    expect_equal(as.double(IC.coxph$ICcumhazard[[1]]), as.double(IC.coxph_d2$ICcumhazard[[1]][order(d2$eventtime),]))
    expect_equal(as.double(IC.cph$ICcumhazard[[1]]), as.double(IClambda_GS[,-1]), tol = 1e-4)
  })
  # }}}
  
  # {{{ before the first event
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

  # {{{ no strata, interactions, continous  
  ICI.coxph <- iidCox(mI.coxph, keep.times = FALSE)
  
  test_that("iid beta - interaction",{
    expect_equal(unname(ICI.coxph$ICbeta),mI.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - interaction",{
    expect_equal(as.double(ICI.coxph$ICcumhazard[[1]]), as.double(IClambdaI_GS[,-1]))
  })
  
    # }}}

  # {{{ no covariate
    IC0.coxph <- iidCox(m0.coxph, keep.times = FALSE) # how to test the result
    # }}}
    
  # {{{ no strata, no interaction, with a categorical variable
  ICCAT.coxph <- iidCox(mCAT.coxph, keep.times = FALSE)
  
  test_that("iid beta - categorical",{
    expect_equal(unname(ICCAT.coxph$ICbeta),mCAT.cox_GS$gamma.iid)
  })
  
  test_that("iid lambda0 - categorical",{
      expect_equal(as.double(ICCAT.coxph$ICcumhazard[[1]]),
                   as.double(IClambdaCAT_GS[,-1]))
  })

    # }}}

  # {{{ no strata but with ties 
  
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

  # {{{ strata, no interaction, continuous
  
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
}

# }}}

# {{{ Cox predict SE

ci2se <- function(x){
    res <- (x$survival.upper-x$survival.lower)/(qnorm(0.975)-qnorm(0.025))
    return(as.double(res))
    }
    
if(require(timereg)){
  
  # {{{ no strata, no interaction, continous    
   
    test_that("predictionsSE",{
        predGS <- predict(m.cox_GS, newdata = d, times = 10)
        predRR1 <- predictCox(m.coxph, newdata = d, times = 10, ci = TRUE)
        predRR2 <- predictCox(m.coxph_d2, newdata = d, times = 10, ci = TRUE)

    
    expect_equal(ci2se(predRR1), as.double(predGS$se.S0))
    expect_equal(ci2se(predRR2), as.double(predGS$se.S0))
    
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e8, ci = TRUE)
    expect_true(all(is.na(predRR1$survival.lower)))
    expect_true(all(is.na(predRR1$survival.upper)))
    
    predRR1 <- predictCox(m.coxph, newdata = d, times = 1e-8, ci = TRUE)
    expect_true(all(predRR1$survival.lower==predRR1$survival.upper))
    
    predRR1 <- predictCox(m.coxph, newdata = d, times = c(1e-8,10,1e8), ci = TRUE)
    expect_true(all(is.na(predRR1$survival.lower[,3])))
    expect_true(all(is.na(predRR1$survival.upper[,3])))
    expect_true(all(predRR1$survival.lower[,1]==predRR1$survival.upper[,1]))
})
  # }}}
  

  # {{{ no strata, interactions, continous  
    test_that("predictionsSE - interaction",{
        predGS <- predict(mI.cox_GS, newdata = d, times = 10)
        predRR1 <- predictCox(mI.coxph, newdata = d, times = 10, ci = TRUE)
        expect_equal(ci2se(predRR1), as.double(predGS$se.S0))
    })
  
    # }}}
    
  # {{{ no strata, no interaction, with a categorical variable
    test_that("predictionsSE - categorical",{
        predGS <- predict(mCAT.cox_GS, newdata = d, times = 10)
        predRR1 <- predictCox(mCAT.coxph, newdata = d, times = 10, ci = TRUE)
        expect_equal(ci2se(predRR1), as.double(predGS$se.S0))
    })

    # }}}

  # {{{ strata, no interaction, continuous

    test_that("predictionsSE - strata",{
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 2)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 2, ci = TRUE)
    
    expect_equal(ci2se(predRR1), as.double(predGS$se.S0))
    
    predGS <- predict(mStrata.cox_GS, newdata = dStrata, times = 1:3)
    predRR1 <- predictCox(mStrata.coxph, newdata = dStrata, times = 1:3, ci = TRUE)
    
    expect_equal(ci2se(predRR1), as.double(predGS$se.S0))
    })
 
    # }}}
}

# }}}



# {{{ Cox confidence band
# package.source("riskRegression", Ccode = TRUE, RorderDescription=FALSE)
# 
#   predGS <- predict(m.cox_GS, newdata = d, times = 10)
#   predRR1 <- predictCox(m.coxph, newdata = d, times = 10, ci = TRUE)
# 
#   bandCox(m.coxph,
#           newdata = d[1:5,],
#           times = 10)
# Edrr2 <- predictCox(m.coxph_d2,
#                     newdata = d[1:5,],
#                     times = 10,
#                     ci = TRUE)


# }}}

##  ### uniform confidence bands, based on resampling  ## {{{
##     if (uniform==1) {
##       mpt <- .C('confBandBasePredict',
##                 delta = as.double(delta), nObs = as.integer(nobs), nt = as.integer(nt),
##                 n = as.integer(n), se = as.double(se), mpt = double(n.sim*nobs),
##                 nSims = as.integer(n.sim), PACKAGE="timereg")$mpt;
  
##       mpt <- matrix(mpt,n.sim,nobs,byrow = TRUE);
##       uband <- apply(mpt,2,percen,per=1-alpha);
##     } else uband<-NULL; 
## ## }}}
