library(timereg)
library(testthat)
library(riskRegression)

n.sim <- 500
# {{{ survival

# {{{ simulation
set.seed(10)
d <- sampleData(1e2, outcome = "survival")
newdata <- d[1:10,]

fit <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X2), data = d, max.timepoint.sim=NULL)
fit.coxph <- coxph(Surv(time, event) ~ X1 + X2, data = d, x = TRUE, y = TRUE)

object.design <- riskRegression:::coxDesign.coxph(fit.coxph)
times <- unique(sort(object.design[object.design$status==1,"stop"]))
# }}}

# {{{ compute quantile for confidence bands
resTimereg <- list()
    for(i in 1:NROW(newdata)){ # i <- 1
        set.seed(10)
        resTimereg[[i]] <- predict.timereg(fit,
                                           newdata = newdata[i,,drop=FALSE],
                                           times = fit$time.sim.resolution,
                                           resample.iid = 1,
                                           n.sim = n.sim)
    }
test_that("computation of the quantile for the confidence band of the cumhazard", {

    set.seed(10)
    ref <- unlist(lapply(resTimereg,"[[", "unif.band"))
    
    pred <- predictCox(fit.coxph,
                       newdata = newdata,
                       times = times,
                       se = TRUE,
                       iid = TRUE,
                       band = TRUE,
                       type = "cumhazard")

    pred.confint <- confint(pred, nsim.band = n.sim, seed = 10)
    expect_equal(pred.confint$cumhazard.quantileBand,ref)

    set.seed(10)
    pred.band2 <- riskRegression:::confBandCox(iid = pred$cumhazard.iid,
                                               se = pred$cumhazard.se,
                                               n.sim = n.sim, 
                                               conf.level = 0.95)

    expect_equal(pred.band2,ref)
})
# }}}

# {{{ display
predRR <- predictCox(fit.coxph,
                     newdata = newdata[1,],
                     times = times,
                     se = TRUE,
                     band = TRUE,
                     type = c("cumhazard","survival")
                     )


plotRR <- autoplot(predRR, type = "survival", band = TRUE, ci = TRUE, plot = FALSE)
dev.new()
plotTR <- plot.predict.timereg(resTimereg[[1]])
dev.new()
plotRR$plot + coord_cartesian(ylim = c(0,1))
# }}}

# {{{ example
test1 <- data.frame(time=c(4,3,1,1,2,2,3), 
                    status=c(1,1,1,0,1,1,0), 
                    x=c(0,2,1,1,1,0,0), 
                    sex=c(0,0,0,0,1,1,1)) 
# Fit a stratified model 
m <- coxph(Surv(time, status) ~ x + strata(sex), 
           data = test1, x = TRUE, y = TRUE) 

set.seed(1)
res <- predictCox(m, newdata = test1, times = 1:4, band = TRUE)
# res$quantile.band
# [1] 2.148709 2.288384 2.327076 2.327076 1.955180 1.951997 1.951997
# }}}

# }}}

# {{{ absolute risk

# {{{ simulation
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")
newdata <- d[1:10,]

fit.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = d)
seqTimes <- fit.CSC$eventTimes

res <- predict(fit.CSC,
               newdata = newdata,
               times = seqTimes-1e-5,
               band = TRUE,
               nsim.band = 500,
               cause = 1)

res <- predict(fit.CSC,
               newdata = newdata[1,,],
               times = seqTimes-1e-5,
               nsim.band = 500,
               band = TRUE, se = TRUE,
               cause = 1)
autoplot(res, band = TRUE, ci = TRUE)


setkey(d,time)
res$times - d$time[d$event>0]
d$event[d$event>0]

# }}}

# }}}


## for debuging

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


