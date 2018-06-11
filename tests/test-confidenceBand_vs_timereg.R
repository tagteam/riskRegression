### test-confidenceBand_vs_timereg.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 18 2017 (09:23) 
## Version: 
## last-updated: jun  6 2018 (18:34) 
##           By: Brice Ozenne
##     Update #: 105
#----------------------------------------------------------------------
## 
### Commentary: 
## Compare the confidence bands returned by predictCox to the one of timereg
##
## TO BE DONE in competing risk setting.
## PB don't know how to use comp.risk to get a CSC.
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(timereg)
library(testthat)
library(riskRegression)

n.sim <- 500
context("[predictCox] confidence band")

## * Survival case

## ** Data
set.seed(10)
dt <- sampleData(1e2, outcome = "survival")
newdata <- dt[1:10,]

dtStrata <- data.frame(time=c(4,3,1,1,2,2,3), 
                       status=c(1,1,1,0,1,1,0), 
                       x=c(0,2,1,1,1,0,0), 
                       sex=c(0,0,0,0,1,1,1)) 

## ** Model
e.timereg <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X2), data = dt, max.timepoint.sim=NULL)
e.coxph <- coxph(Surv(time, event) ~ X1 + X2, data = dt, x = TRUE, y = TRUE)

vec.times <- e.timereg$time.sim.resolution

## ** Compute quantile for confidence bands
resTimereg <- list()
for(i in 1:NROW(newdata)){ # i <- 1
    set.seed(10)
    resTimereg[[i]] <- predict.timereg(e.timereg,
                                       newdata = newdata[i,,drop=FALSE],
                                       times = vec.times,
                                       resample.iid = 1,
                                       n.sim = n.sim)
}

## ** Tests
test_that("[predictCox] Quantile for the confidence band of the cumhazard", {

    predRR <- predictCox(e.coxph,
                         newdata = newdata,
                         times = vec.times,
                         se = TRUE,
                         iid = TRUE,
                         band = TRUE,
                         type = "cumhazard")

    ## compatibility with timereg
    ref <- unlist(lapply(resTimereg,"[[", "unif.band"))

    set.seed(10)
    pred.band2 <- riskRegression:::confBandCox(iid = predRR$cumhazard.iid,
                                               se = predRR$cumhazard.se,
                                               n.sim = n.sim, 
                                               conf.level = 0.95)

    expect_equal(pred.band2,ref)


    ## note confint is removing the first column since the standard error is 0
    set.seed(10)
    pred.band2.no0 <- riskRegression:::confBandCox(iid = predRR$cumhazard.iid[,-1,,drop=FALSE],
                                                   se = predRR$cumhazard.se[,-1],
                                                   n.sim = n.sim, 
                                                   conf.level = 0.95)

    ## should not set transform to NA because at time 0 se=0 so the log-transform fails
    pred.confint <- confint(predRR, nsim.band = n.sim, seed = 10,
                            cumhazard.transform = "none")
    expect_equal(pred.confint$cumhazard.quantileBand, pred.band2.no0)
    expect_equal(pred.confint$cumhazard.quantileBand, ref)

})

## ** Display
predRR <- predictCox(e.coxph,
                     newdata = newdata[1],
                     times = vec.times,
                     se = TRUE,
                     band = TRUE,
                     type = c("cumhazard","survival")
                     )


plotRR <- autoplot(predRR, type = "survival", band = TRUE, ci = TRUE, plot = FALSE)
dev.new()
plotTR <- plot.predict.timereg(resTimereg[[1]])
dev.new()
plotRR$plot + coord_cartesian(ylim = c(0,1))
graphics.off()

## ** With strata                                        
## Fit a stratified model 
eS.coxph <- coxph(Surv(time, status) ~ x + strata(sex), 
                  data = dtStrata, x = TRUE, y = TRUE) 

eS.pred  <- predictCox(eS.coxph, newdata = dtStrata, times = 1:4, band = TRUE)
eS.confint <- confint(eS.pred, seed = 10)
eS.confint$survival.quantileBand


## * Competing risk setting

## ** Data
set.seed(10)
dt <- sampleData(1e2, outcome = "competing.risks")
newdata <- dt[1:10,]

## ** Model
e.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = dt)
vec.times <- e.CSC$eventTimes

e.timereg <- comp.risk(Event(time, event) ~ const(X1) + const(X2),
                       model = "rcif",
                       data = dt, cause = 1)
coef(e.timereg)
coef(e.CSC)

## ** Compute confidence bands
if(FALSE){
    resTimereg <- predict.timereg(e.timereg,
                                  newdata = newdata[1,,drop=FALSE],
                                  times = vec.times,
                                  resample.iid = 1,
                                  n.sim = n.sim)
    resTimereg$P1
}
predRR <- predict(e.CSC,
                  newdata = newdata,
                  times = vec.times-1e-5,
                  se = TRUE,
                  band = TRUE,
                  nsim.band = 500,
                  cause = 1)

predRR$absRisk[1,]

## ** Display
autoplot(predRR, band = TRUE, ci = TRUE)
graphics.off()

## * for debuging

## cumHazard.coxph <- predictCox(fit.coxph)$cumhazard
## iid.coxph <- iidCox(fit.coxph)
## iid.lambda <- iid.coxph$ICcumhazard[[1]][1,]

## X.design <- model.matrix(formula(fit.coxph),newdata[i,,drop=FALSE])[,-1]
## eLP <- exp(X.design %*% cbind(coef(fit.coxph)))
## Xiid.beta <- X.design %*% iid.coxph$IFbeta[1,]

## term1 <- as.numeric(eLP * iid.lambda)
## term2 <- as.numeric(eLP * cumHazard.coxph * Xiid.beta)

## ls.args <- list(delta = res$cumhazard.iid[i,,],
##                 nObs = 1,
##                 nt = length(times),
##                 n = NROW(d),
##                 mpt = n.sim*1,
##                 nSims = n.sim)
## ls.args$se <-  apply(ls.args$delta^2, 1, sum)^0.5

## mpt <- .C("confBandBasePredict", delta = as.double(delta), 
##           nObs = as.integer(nobs), nt = as.integer(nt), n = as.integer(n), 
##           se = as.double(se), mpt = double(n.sim * nobs), nSims = as.integer(n.sim), 
##           PACKAGE = "timereg")$mpt


#----------------------------------------------------------------------
### test-confidenceBand_vs_timereg.R ends here
