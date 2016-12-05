library(prodlim)
library(riskRegression)
library(testthat)
library(rms)
library(survival)

####
context("Cox prediction - iid")

if(require(timereg)){
  
  set.seed(10)
  d <- sampleData(5e1, outcome = "survival")[,.(eventtime,event,X1,X2,X6)]
  d[ , X16 := X1*X6]
  d2 <- copy(d)
  setkey(d,eventtime)
  
  #### non stratified Cox model
  mGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  IClambda_timereg <- t(as.data.table(mGS.cox$B.iid))
  
  m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE)
  IC.cox <- iidCox(m.cox, keep.times = FALSE)
  
  m.cox2 <- coxph(Surv(eventtime, event) ~ X1+X6, data = d2, y = TRUE)
  IC.cox2 <- iidCox(m.cox2, keep.times = FALSE) #### bug here when n = 50
  
  m.cox3 <- cph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE)
  IC.cox3 <- iidCox(m.cox3, keep.times = FALSE)
  
  test_that("iid beta",{
    expect_equal(IC.cox$ICbeta,mGS.cox$gamma.iid)
    expect_equal(IC.cox$ICbeta[order(d2$eventtime),],IC.cox2$ICbeta)
    expect_equal(IC.cox3$ICbeta,mGS.cox$gamma.iid, tol = 1e-2)
  })
  
  test_that("iid lambda0",{
    expect_equal(as.double(IC.cox$ICLambda0), as.double(IClambda_timereg[,-1]))
    expect_equal(as.double(IC.cox$ICLambda0[order(d2$eventtime),]), as.double(IC.cox2$ICLambda0))
    expect_equal(as.double(IC.cox3$ICLambda0), as.double(IClambda_timereg[,-1]), tol = 1e-4)
  })
  
  #### non stratified Cox model with interactions
  mIGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1) + prop(X6) + prop(X1*X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  ICIlambda_timereg <- t(as.data.table(mIGS.cox$B.iid))
  
  mI.cox <- coxph(Surv(eventtime, event) ~ X1*X6, data = d, y = TRUE)
  IC.Icox <- iidCox(mI.cox, keep.times = FALSE)
  
  test_that("iid beta - interaction",{
    expect_equal(IC.Icox$ICbeta,mIGS.cox$gamma.iid)
  })
  
  test_that("iid lambda0 - interaction",{
    expect_equal(as.double(IC.Icox$ICLambda0), as.double(ICIlambda_timereg[,-1]))
  })
  
  #### Cox model with no covariate
  # mIGS.cox <- cox.aalen(Surv(eventtime, event) ~ 1, data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  # ICIlambda_timereg <- t(as.data.table(mIGS.cox$B.iid))
  
  mI.cox <- coxph(Surv(eventtime, event) ~ 1, data = d, y = TRUE)
  IC.Icox <- iidCox(mI.cox, keep.times = FALSE) # how to test the result
  
  #### Cox model with categorical variable
  d$Xcat2 <- as.factor(paste0(d$X1,d$X2))
  
  m.timereg <- cox.aalen(Surv(eventtime, event) ~ prop(Xcat2) + prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  IClambda.timereg <- t(as.data.table(m.timereg$B.iid))
  
  m.RR <- coxph(Surv(eventtime, event) ~ Xcat2 + X6, data = d, y = TRUE)
  IC.RR <- iidCox(m.RR, keep.times = FALSE)
  
  test_that("iid beta - categorical",{
    expect_equal(IC.RR$ICbeta,m.timereg$gamma.iid)
  })
  
  test_that("iid lambda0 - categorical",{
    expect_equal(as.double(IC.RR$ICLambda0), as.double(IClambda.timereg[,-1]))
  })
  
  #### stratified Cox model
  # dStrata <- rbind(cbind(d, St= 1),
  #                  cbind(d, St= 2))
  # mSGS.cox <- cox.aalen(Surv(eventtime, event) ~ strata(St) + prop(X1) + prop(X6), data = dStrata, resample.iid = TRUE, max.timepoint.sim=NULL)
  # ICSlambda_timereg <- t(as.data.table(mSGS.cox$B.iid))
  # 
  # setkeyv(dStrata, c("St", "eventtime"))
  # mS.cox <- coxph(Surv(eventtime, event) ~ strata(St) + X1 + X6, data = dStrata, y = TRUE)
  # IC.Scox <- iidCox(mS.cox, keep.times = FALSE)
  # 
  # crossprod(mSGS.cox$gamma.iid) / crossprod(IC.Scox$ICbeta)
  # cbind(mSGS.cox$gamma.iid,IC.Scox$ICbeta)
}

####
context("Cox prediction - standard error")


set.seed(10)
d <- sampleData(1e2, outcome = "survival")
nd <- sampleData(100, outcome = "survival")
d$time <- round(d$time,1)
d <- d[order(d$time),]

for(ties in c("breslow","efron")){ # ties <- "breslow"
  test_that(paste("predictCox(empty) - valide se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ 1,data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ 1,data=d, method=ties, y = TRUE)

    # res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = FALSE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = FALSE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
  })
  
  test_that(paste("predictCox(univariate) - valid se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ X1,data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ X1,data=d, method=ties, y = TRUE, surv = TRUE)
    
    res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
    
    expect_equal(res_surv$fit,diag(resCoxph$cumHazard))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
    # expect_equal(resCph, resCoxph, tolerance = 1e-4, scale = 1) # expected difference coef(fit_cph)-coef(fit_coxph)
    # pec:::predictSurvProb(fit_coxph, newdata = d, times = d$time) - pec:::predictSurvProb(fit_cph, newdata = d, times = d$time)
  })
  
  test_that(paste("predictCox(multiple) - valid se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ X1 + X2 + X6,data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ X1 + X2 + X6,data=d, method=ties, y = TRUE)
    
    res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
    res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
    
    expect_equal(res_surv$fit,diag(resCoxph$cumHazard))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))  # expected difference coef(fit_cph)-coef(fit_coxph)
    # expect_equal(resCph, resCoxph, tolerance = 1e-5, scale = 1)
  })
}

#### strata
for(ties in c("breslow","efron")){ #
    test_that(paste("predictCox (strata, empty) - valide se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ strata(X1),data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ strat(X1),data=d, method=ties, y = TRUE)

    # res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = FALSE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = FALSE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
    # res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
  })
  
  test_that(paste("predictCox (strata, univariate) - valide se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ strata(X1) + X2,data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ strat(X1) + X2,data=d, method=ties, y = TRUE)
    
    res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
    res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
    
    expect_equal(res_surv$fit,diag(t(resCoxph$cumHazard)))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
    # expect_equal(resCph,resCoxph, tolerance = 1e-5, scale = 1)
  })
  
  test_that(paste("predictCox (strata, multiple) - valide se cumHazard",ties),{
    fit_coxph <- coxph(Surv(time,event) ~ strata(X1) + X2 + X6 ,data=d, ties=ties)
    fit_cph <- cph(Surv(time,event) ~ strat(X1) + X2 + X6 ,data=d, method=ties, y = TRUE)
    
    res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
    resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
    resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
    res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
    
    expect_equal(res_surv$fit,diag(t(resCoxph$cumHazard)))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
    # expect_equal(resCph,resCoxph, tolerance = 1e-5, scale = 1)
  })
}

#### bootstrap version

test <- FALSE
if(test){
  
  n <- 1e3
  df <- sampleData(n, outcome =  "survival")
  df.test <- df[1:10,]
  
  fit_coxph <- coxph(Surv(time,event) ~ X1 + X2 + X3 + X4 + X5, data=df, ties="efron", y = TRUE)
  resPredict <- predictCox(fit_coxph, newdata = df.test, times = df.test$time, se = TRUE)
  resSurv <- survival:::predict.coxph(fit_coxph, newdata = df.test, type = "expected", se.fit = TRUE)
  expect_equal(diag(resPredict$cumHazard), resSurv$fit)
  expect_equal(diag(resPredict$cumHazard.se), resSurv$se.fit)
  
  n.bootstrap <- 1e3
  array.boot <- list(hazard = array(NA, dim = c(nrow(df.test), nrow(df.test), n.bootstrap)),
                     cumHazard = array(NA, dim = c(nrow(df.test), nrow(df.test), n.bootstrap)),
                     survival = array(NA, dim = c(nrow(df.test), nrow(df.test), n.bootstrap))
  )
  
  pb <- utils::txtProgressBar(min = 0, max = n.bootstrap, initial = 0)
  for(iterB in 1:n.bootstrap){
    
    data_boot <- df[sample.int(nrow(df), nrow(df), replace = TRUE),]
    fit_boot <- coxph(Surv(time,event) ~ X1 + X2 + X3 + X4 + X5, data=data_boot, 
                      ties="efron", y = TRUE)
    res <- predictCox(fit_boot, newdata = df.test, times = df.test$time, se = FALSE)
    
    array.boot$hazard[,,iterB] <- res$hazard
    array.boot$cumHazard[,,iterB] <- res$cumHazard
    array.boot$survival[,,iterB] <- res$survival
    setTxtProgressBar(pb, iterB)
  }
  close(pb)
  
  ls.seBOOT <- list(hazard.se = apply(array.boot$hazard, 1:2, sd),
                    cumHazard.se = apply(array.boot$cumHazard, 1:2, sd),
                    survival.se = apply(array.boot$survival, 1:2, sd)
  )
  
  ## NOT GOOD!!!
  (ls.seBOOT$cumHazard.se-resPredict$hazard.se)/ls.seBOOT$hazard.se
  (ls.seBOOT$cumHazard.se-resPredict$cumHazard.se)/ls.seBOOT$cumHazard.se
  (ls.seBOOT$cumHazard.se-resPredict$survival.se)/ls.seBOOT$survival.se
  
}

