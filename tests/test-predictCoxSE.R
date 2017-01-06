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
  
  m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
  IC.cox <- iidCox(m.cox, keep.times = FALSE)
  
  m.cox2 <- coxph(Surv(eventtime, event) ~ X1+X6, data = d2, y = TRUE, x = TRUE)
  IC.cox2 <- iidCox(m.cox2, keep.times = FALSE)
  
  m.cox3 <- cph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
  IC.cox3 <- iidCox(m.cox3, keep.times = FALSE)
  
  test_that("iid beta",{
    expect_equal(unname(IC.cox$ICbeta),mGS.cox$gamma.iid)
    expect_equal(unname(IC.cox$ICbeta),unname(IC.cox2$ICbeta[order(d2$eventtime),]))
    expect_equal(unname(IC.cox3$ICbeta),mGS.cox$gamma.iid, tol = 1e-2)
  })
  
  test_that("iid lambda0",{
    expect_equal(as.double(IC.cox$ICcumHazard[[1]]), as.double(IClambda_timereg[,-1]))
    expect_equal(as.double(IC.cox$ICcumHazard[[1]]), as.double(IC.cox2$ICcumHazard[[1]][order(d2$eventtime),]))
    expect_equal(as.double(IC.cox3$ICcumHazard[[1]]), as.double(IClambda_timereg[,-1]), tol = 1e-4)
  })
  
  test_that("predictionsSE",{
    predGS <- predict(mGS.cox, newdata = d, times = 10)
    predRR1 <- predictCox(m.cox, newdata = d, times = 10, se = TRUE)
    predRR2 <- predictCox(m.cox, newdata = d, times = 10, se = TRUE, iid = iidCox(m.cox))
    expect_equal(predRR1$survival.se, predGS$se.S0)
    expect_equal(predRR2$survival.se, predGS$se.S0)
    
    predRR1 <- predictCox(m.cox, newdata = d, times = 1e8, se = TRUE)
    expect_true(all(is.na(predRR1$survival.se)))
    
    predRR1 <- predictCox(m.cox, newdata = d, times = 1e-8, se = TRUE)
    expect_true(all(predRR1$survival.se==0))
    
    predRR1 <- predictCox(m.cox, newdata = d, times = c(1e-8,10,1e8), se = TRUE)
    expect_true(all(is.na(predRR1$survival.se[,3])))
    expect_true(all(predRR1$survival.se[,1]==0))
  })
  
  #### before the first event
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
    expect_equal(as.double(iid1$ICcumHazard[[1]]), as.double(ICGS_lambda0[,-1]))
  })
  
  #### non stratified Cox model with interactions
  mIGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1) + prop(X6) + prop(X1*X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  ICIlambda_timereg <- t(as.data.table(mIGS.cox$B.iid))
  
  mI.cox <- coxph(Surv(eventtime, event) ~ X1*X6, data = d, y = TRUE, x = TRUE)
  IC.Icox <- iidCox(mI.cox, keep.times = FALSE)
  
  test_that("iid beta - interaction",{
    expect_equal(unname(IC.Icox$ICbeta),mIGS.cox$gamma.iid)
  })
  
  test_that("iid lambda0 - interaction",{
    expect_equal(as.double(IC.Icox$ICcumHazard[[1]]), as.double(ICIlambda_timereg[,-1]))
  })
  
  test_that("predictionsSE - interaction",{
    predGS <- predict(mIGS.cox, newdata = d, times = 10)
    predRR1 <- predictCox(mI.cox, newdata = d, times = 10, se = TRUE)
    expect_equal(predRR1$survival.se, predGS$se.S0)
  })
  
  #### Cox model with no covariate
  # mIGS.cox <- cox.aalen(Surv(eventtime, event) ~ 1, data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  # ICIlambda_timereg <- t(as.data.table(mIGS.cox$B.iid))
  
  mI.cox <- coxph(Surv(eventtime, event) ~ 1, data = d, y = TRUE, x = TRUE)
  IC.Icox <- iidCox(mI.cox, keep.times = FALSE) # how to test the result
  
  #### Cox model with categorical variable
  d$Xcat2 <- as.factor(paste0(d$X1,d$X2))
  
  m.timereg <- cox.aalen(Surv(eventtime, event) ~ prop(Xcat2) + prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
  IClambda.timereg <- t(as.data.table(m.timereg$B.iid))
  
  m.RR <- coxph(Surv(eventtime, event) ~ Xcat2 + X6, data = d, y = TRUE, x = TRUE)
  IC.RR <- iidCox(m.RR, keep.times = FALSE)
  
  test_that("iid beta - categorical",{
    expect_equal(unname(IC.RR$ICbeta),m.timereg$gamma.iid)
  })
  
  test_that("iid lambda0 - categorical",{
    expect_equal(as.double(IC.RR$ICcumHazard[[1]]), as.double(IClambda.timereg[,-1]))
  })
  
  test_that("predictionsSE - interaction",{
    predGS <- predict(m.timereg, newdata = d, times = 10)
    predRR1 <- predictCox(m.RR, newdata = d, times = 10, se = TRUE)
    expect_equal(predRR1$survival.se, predGS$se.S0)
  })
  
  #### Ties
  d3 <- copy(d)[1:10,]
  d3[, event := 1]
  d3[7:8, X1 := 1]
  d3[7:8, X6 := 1]
  setkey(d3, eventtime)
  d3[, eventtimeTies := eventtime]
  d3[7:8, eventtimeTies := eventtime[1]]
  # d3[, eventtime := round(eventtime, 0)]
  
  mGS.cox0 <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d3, resample.iid = TRUE, max.timepoint.sim=NULL)
  mGS.cox <- cox.aalen(Surv(eventtimeTies, event) ~ prop(X1)+prop(X6), data = d3, resample.iid = TRUE, max.timepoint.sim=NULL)
  IClambda_timereg <- t(as.data.table(mGS.cox$B.iid))
  
  m.cox0 <- coxph(Surv(eventtime, event) ~ X1+X6, data = d3, y = TRUE, x = TRUE)
  m.cox <- coxph(Surv(eventtimeTies, event) ~ X1+X6, data = d3, y = TRUE, x = TRUE)
  IC.cox0 <- iidCox(m.cox0, keep.times = FALSE)
  IC.cox <- iidCox(m.cox, keep.times = FALSE)
  
  # test_that("iid beta - categorical",{
  #   expect_equal(unname(IC.cox$ICbeta),mGS.cox$gamma.iid)
  # })
  data.frame(RR0 = mGS.cox0$gamma.iid, GS0 = IC.cox0$ICbeta, RR = IC.cox$ICbeta, GS = mGS.cox$gamma.iid)
  
  # test_that("iid lambda0 - categorical",{
  #   expect_equal(as.double(IC.RR$ICcumHazard[[1]]), as.double(IClambda.timereg[,-1]))
  # })
  
  #### stratified Cox model
  # dStrata <- rbind(cbind(d[1:10], St= 1),
  #                  cbind(d[1:10], St= 2))
  dStrata <- d
  dStrata$St <- rbinom(n = NROW(d), size = 2, prob = c(1/3,1/2))
  dStrata2 <- copy(dStrata)
  setkeyv(dStrata, c("St", "eventtime"))
  
  mGSS.cox <- cox.aalen(Surv(eventtime, event) ~ strata(St)-1 + prop(X1) + prop(X6), data = dStrata, 
                        resample.iid = TRUE, max.timepoint.sim=NULL)
  
  mS.cox <- coxph(Surv(eventtime, event) ~ strata(St) + X1 + X6, data = dStrata, y = TRUE, x = TRUE)
  IC.Scox <- iidCox(mS.cox)
  
  test_that("iid beta - strata",{
    expect_equal(unname(IC.Scox$ICbeta),mGSS.cox$gamma.iid)
  })
  
  test_that("iid lambda0 - strata",{
    
    for(iStrata in 1:length(unique(dStrata$St))){
      IC.GS <- do.call(rbind,
                       lapply(mGSS.cox$B.iid,function(x){x[,iStrata]})
      )
      colnames(IC.GS) <- mGSS.cox$time.sim.resolution
      
      checkTimes <- intersect(mGSS.cox$time.sim.resolution,IC.Scox$time[[iStrata]])
      
      
      diff <- IC.Scox$ICcumHazard[[iStrata]][,which(IC.Scox$time[[iStrata]] %in% checkTimes),drop = FALSE]-IC.GS[,which(mGSS.cox$time.sim.resolution %in% checkTimes)]
      expect_true(all(abs(na.omit(as.double(diff)))<1e-10))
    }
    
  })
  
  test_that("predictionsSE - strata",{
    predGS <- predict(mGSS.cox, newdata = dStrata, times = 2)
    predRR1 <- predictCox(mS.cox, newdata = dStrata, times = 2, se = TRUE)
    
    expect_equal(predRR1$survival.se, predGS$se.S0)
    
    predGS <- predict(mGSS.cox, newdata = dStrata, times = 1:3)
    predRR1 <- predictCox(mS.cox, newdata = dStrata, times = 1:3, se = TRUE)
    
    expect_equal(predRR1$survival.se, predGS$se.S0)
  })
  
}


#### bootstrap version
test <- FALSE
if(test){
  
  set.seed(10)
  d <- SimCompRisk(1e3)
  d <- d[order(d$time),c("time", "event", "X1", "X2", "cause")]
  ttt <- sample(x = unique(sort(d$time)), size = 3)
  d <- d[order(d$time), ]
  
  library(boot)
  
  #### survival
  predCox <- function(d,i){
    coxph.fit <- coxph(Surv(time,event==1)~ X1+X2,data=d[i,], method = "breslow", x = TRUE, y = TRUE)
    res <- predictCox(coxph.fit, newdata = d[1:2,,drop=FALSE], times = seq(2,5,1), se = FALSE, type = "survival")
    return(res$survival)
  }
  predCox(d, 1:NROW(d))
  res.boot <- boot(d, predCox, R = 1000, stype = "i", ncpus = 4)
  
  coxph.fit <- coxph(Surv(time,event==1)~ X1+X2,data=d, method = "breslow", x = TRUE, y = TRUE)
  res.IF <- predictCox(coxph.fit, newdata = d[1:2,,drop=FALSE], times = seq(2,5,1), se = TRUE, type = "survival")
  
  res.boot$t0-res.IF$survival
  apply(res.boot$t,2,sd) - as.double(res.IF$survival.se)
  (apply(res.boot$t,2,sd) - as.double(res.IF$survival.se))/as.double(res.IF$survival.se)
  
  #### absolute risk
  predCSC <- function(d,i){
    CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d[i,], method = "breslow", iid = FALSE)
    res <- predict(CSC.fit, newdata = d[1:2,,drop=FALSE], cause = 1, times = seq(2,5,1), se = FALSE)
    return(res)
  }
  
  system.time(
    res.boot <- boot(d, predCSC, R = 500, stype = "i", ncpus = 4)
  )
  
  CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", iid = TRUE)
  res.IF <- predict(CSC.fit, newdata = d[1:2,,drop=FALSE], cause = 1, times = seq(2,5,1), se = TRUE)
  
  res.boot$t0-res.IF
  apply(res.boot$t,2,sd) - as.double(attr(res.IF,"se"))
  (apply(res.boot$t,2,sd) - as.double(attr(res.IF,"se")))/as.double(attr(res.IF,"se"))
  
  
}





