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

object.design <- riskRegression:::CoxDesign.coxph(fit.coxph)
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
    pred <- predictCox(fit.coxph,
                       newdata = newdata,
                       times = times,
                       se = TRUE,
                       iid = TRUE,
                       band = TRUE,
                       type = "cumhazard")
    set.seed(10)
    resRR <- riskRegression:::confBandCox(iid = pred$cumhazard.iid,
                                          se = pred$cumhazard.se,
                                          times = times,
                                          n.sim = n.sim, conf.level = 0.95)


    
    ref <- unlist(lapply(resTimereg,"[[", "unif.band"))

    set.seed(10)
    predRR <- predictCox(fit.coxph,
                         newdata = newdata,
                         times = times,
                         se = TRUE,
                         band = TRUE,
                         nSim.band = n.sim,
                         type = c("cumhazard")
                         )
    expect_equal(predRR$quantile.band,ref)
    expect_equal(resRR,ref)
    #expect_equal(resRR,c(2.850666, 2.562949, 2.992631, 2.992631, 2.850666, 2.992631, 2.850666, 2.992631, 2.992631, 2.562949), tol = 1e-6)
    #expect_equal(predRR$quantile.band,c(2.850666, 2.562949, 2.992631, 2.992631, 2.850666, 2.992631, 2.850666, 2.992631, 2.992631, 2.562949), tol = 1e-6)
})
# }}}

# {{{ display
## predRR <- predictCox(fit.coxph,
                     ## newdata = newdata[1,],
                     ## times = times,
                     ## se = TRUE,
                     ## band = TRUE,
                     ## nSim.band = 500,
                     ## type = c("cumhazard","survival")
                     ## )

## dev.new()
## plotRR <- autoplot(predRR, type = "survival", band = TRUE, ci = TRUE, plot = FALSE)

## dev.new()
## plotTR <- plot.predict.timereg(resTimereg[[1]])
## dev.new()
## plotRR$plot + coord_cartesian(ylim = c(0,1))

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
               nSim.band = 500,
               cause = 1)

res <- predict(fit.CSC,
               newdata = newdata[1,,],
               times = seqTimes-1e-5,
               nSim.band = 500,
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
