
library(prodlim)
library(riskRegression)
library(testthat)
library(rms)
library(survival)
context("Cox prediction - standard error")


set.seed(10)
d <- SimSurv(1e2)
nd <- SimSurv(10)
d$time <- round(d$time,1)
d <- d[order(d$time),][1:10,]

for(ties in c("efron")){ #"breslow"
  fit_coxph <- coxph(Surv(time,status) ~ X1 + X2,data=d, ties=ties)
  fit_cph <- cph(Surv(time,status) ~ X1 + X2,data=d, method=ties, y = TRUE)
  # table(duplicated(d$time))
  
  res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
  resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
  
  resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
  res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
  
  test_that(paste("predictCox - valide se cumHazard",ties),{
    expect_equal(res_surv$fit,diag(resCoxph$cumHazard))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
    expect_equal(resCph, resCoxph, tolerance = 1e-5, scale = 1)
  })
}

#### strata
for(ties in c("efron")){ #"breslow"
  fit_coxph <- coxph(Surv(time,status) ~ strata(X1) + X2,data=d, ties=ties)
  fit_cph <- cph(Surv(time,status) ~ strat(X1) + X2,data=d, method=ties, y = TRUE)
  # table(duplicated(d$time))
  
  res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
  resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
  resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
  res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
  
  test_that(paste("predictCox (strata) - valide se cumHazard",ties),{
    expect_equal(res_surv$fit,diag(t(resCoxph$cumHazard)))
    expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
    expect_equal(resCph,resCoxph, tolerance = 1e-5, scale = 1)
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

