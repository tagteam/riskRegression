library(data.table)
library(rbenchmark)
library(prodlim)
library(rms)
library(testthat)
library(survival)
#library(riskRegression)

#### 1- Without strata ####
set.seed(10)
n <- 300
df.NS <- SimCompRisk(n)
df.NS$time <- round(df.NS$time,2)
# table(duplicated(df.NS$time))

seqTime <- c(unique(sort(df.NS$time)), max(df.NS$time) + 1)

for(model in c("coxph","cph")){
  for(method.ties in c("breslow","efron")){
    
    CSC.NS <- CSC(Hist(time,event) ~ X1*X2,data = df.NS, ties = method.ties, fitter = model, x = TRUE, y = TRUE )
    
      for(cause in 1:2){ 
        Event0.NS <- pec:::predictEventProb.CauseSpecificCox(CSC.NS, newdata = df.NS, times = seqTime, cause = cause)
        EventTest.NS <- predictEvent(CSC.NS, newdata = df.NS, times = seqTime, cause = cause)
        
        test_that(paste0("pec vs predictEventProbRR - no strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
          expect_equal(as.numeric(Event0.NS),as.numeric(unlist(EventTest.NS)), tolerance = 1e-8)
        })
      }
  }
}

#### 2- With strata ####
set.seed(10)
n <- 300
df.S <- SimCompRisk(n)
df.S$time <- round(df.S$time,2)
df.S$X3 <- rbinom(n, size = 4, prob = rep(0.25,4))

seqTime <- c(unique(sort(df.S$time)), max(df.S$time) + 1)

for(model in c("coxph","cph")){
  for(method.ties in c("breslow","efron")){
    
    if(model == "coxph"){
      CSC.S <- CSC(Hist(time,event) ~ strata(X1) + strata(X3) + X2, data = df.S, ties = method.ties, fitter = "coxph")
    }else if(model == "cph"){
      CSC.S <- CSC(Hist(time,event) ~ strat(X1) + strat(X3) + X2, data = df.S, ties = method.ties, fitter = "cph", x = TRUE, y = TRUE )
    }
    
    for(cause in 1:2){ 
      Event0.S <- pec:::predictEventProb.CauseSpecificCox(CSC.S, newdata = df.S, times = seqTime, cause = cause)
      EventTest.S <- predictEvent(CSC.S, newdata = df.S, times = seqTime, cause = cause)
      
      test_that(paste0("pec vs predictEventProbRR - strata, method.ties/method.baseHaz/model/cause = ",method.ties,"/",model,"/",cause,""),{
        expect_equal(as.numeric(Event0.S),as.numeric(unlist(EventTest.S)), tolerance = 1e-8)
      })
    }
  }
}
