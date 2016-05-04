library(data.table)
library(rbenchmark)
library(prodlim)
library(rms)
library(testthat)
library(survival)
#library(riskRegression)


#### 1- Without strata ####
set.seed(10)
n <- 1000
df.NS <- SimSurv(n)
df.NS$time <- round(df.NS$time,1)

seq_x <- c(seq(0, max(df.NS$time), length.out = 10), max(df.NS$time) + 1)

for(model in c("coxph","cph")){
  for(method.ties in c("breslow","efron")){
   
     if(model == "coxph"){
      coxph.NS <- coxph(formula = Surv(time,status) ~ X1 * X2, data = df.NS, ties = method.ties)
      Surv0.NS <- pec:::predictSurvProb.coxph(coxph.NS, newdata = df.NS, times = seq_x)
    }else{
      coxph.NS <- cph(formula = Surv(time,status) ~ X1 * X2, data = df.NS, ties = method.ties, surv = TRUE, x = TRUE, y = TRUE)
      Surv0.NS <- pec:::predictSurvProb.cph(coxph.NS, newdata = df.NS, times = seq_x)
    }
    
      SurvTest.NS <- predictSurv(coxph.NS, seq_x, newdata = df.NS)
      HazTest.NS <- predictSurv(coxph.NS, seq_x, newdata = df.NS, type = "hazard")
      CumHazTest.NS <- predictSurv(coxph.NS, seq_x, newdata = df.NS, type = "cumHazard")
      
      test_that(paste0("pec vs predictSurvProbRR - no strata, method.ties/method.baseHaz/model = ",method.ties,"/",model,""),{
        expect_equal(object = as.vector(Surv0.NS), expected = as.numeric(unlist(SurvTest.NS)), tolerance=1e-8)
        expect_equal(object = as.vector(Surv0.NS), expected = as.numeric(unlist(exp(-CumHazTest.NS))), tolerance=1e-8)
      })

  }
}

#### 2- With strata ####
set.seed(10)
n <- 1000

df.S <- SimSurv(n)
df.S$time <- round(df.S$time,1)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))
# table(duplicated(df.S$time))

seq_x <- c(seq(0, max(df.S$time), length.out = 10), max(df.S$time) + 1)

for(model in c("coxph","cph")){
  for(method.ties in c("breslow","efron")){
    
    if(model == "coxph"){
      coxph.S <- coxph(formula = Surv(time,status) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties)
      Surv0.S <- pec:::predictSurvProb.coxph(coxph.S, newdata = df.S, times = seq_x)
    }else{
      coxph.S <- cph(formula = Surv(time,status) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, surv = TRUE, x = TRUE, y = TRUE)
      Surv0.S <- pec:::predictSurvProb.cph(coxph.S, newdata = df.S, times = seq_x)
    }
    
    for(method.baseHaz in c("dt","cpp")){
      
      SurvTest.S <- predictSurv(coxph.S, seq_x, newdata = df.S)
      HazTest.S <- predictSurv(coxph.S, seq_x, newdata = df.S, type = "hazard")
      CumHazTest.S <- predictSurv(coxph.S, seq_x, newdata = df.S, type = "cumHazard")
      
      test_that(paste0("pec vs predictSurvProbRR - strata, method.ties/method.baseHaz/model = ",method.ties,"/",model,""),{
        expect_equal(object = as.vector(Surv0.S), expected = as.numeric(unlist(SurvTest.S)), tolerance = 1e-8)
        expect_equal(object = as.vector(Surv0.S), expected = as.numeric(unlist(exp(-CumHazTest.S))), tolerance = 1e-8)
      })
      
    }
  }
}