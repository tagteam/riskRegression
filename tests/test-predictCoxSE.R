
library(testthat)
library(riskRegression)
library(testthat)
library(rms)
library(survival)
context("Cox prediction - standard error")


library(rms)
set.seed(10)
d <- SimSurv(1e2)
nd <- SimSurv(10)
d$time <- round(d$time,1)
d <- d[order(d$time),][1:10,]
fit_coxph <- coxph(Surv(time,status) ~ X1 + X2,data=d, ties="efron")
fit_cph <- cph(Surv(time,status) ~ X1 + X2,data=d, ties="efron", y = TRUE)
# table(duplicated(d$time))

res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
expect_equal(res_surv$fit,diag(t(resCoxph$cumHazard)))
expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
expect_equal(resCph,resCoxph, tolerance = 1e-5, scale = 1)



fit_coxph <- coxph(Surv(time,status) ~ strata(X1) + X2,data=d, ties="efron")
fit_cph <- cph(Surv(time,status) ~ strat(X1) + X2,data=d, ties="efron", y = TRUE)
# table(duplicated(d$time))

res_surv <- survival:::predict.coxph(fit_coxph, newdata = d, type="expected", se.fit = TRUE)
resCoxph <- predictCox(fit_coxph, newdata = d, times = d$time,  se = TRUE)
resCph <- predictCox(fit_cph, newdata = d, times = d$time,  se = TRUE)
res_pec <- pec::predictSurvProb(fit_coxph, newdata = d, times = d$time)
expect_equal(res_surv$fit,diag(t(resCoxph$cumHazard)))
expect_equal(res_surv$se.fit,diag(resCoxph$cumHazard.se))
expect_equal(resCph,resCoxph, tolerance = 1e-5, scale = 1)
