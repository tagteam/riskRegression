# library(butils.base)

library(testthat)
library(riskRegression)
library(survival)
library(rms)

set.seed(10)
n <- 20
d <- SimCompRisk(n)
d <- d[order(d$time),c("time", "event", "X1", "X2", "cause")]
timeFirstEvent <- min(d[d$event==1,]$time)
times <-  timeFirstEvent+1:5
type <- c("survival","cumhazard")
newdata <- d[1:4,]

range(d$time)

coxph.fit <- coxph(Surv(time,event==1)~ X1+X2,data=d, method = "breslow", x = TRUE, y = TRUE)

resSE <- predictCox(coxph.fit, newdata = newdata, times = times, se = TRUE, type = type)
resIID <- predictCox(coxph.fit, newdata = newdata, times = times, iid = TRUE, type = type)

# {{{ "Export IID for survival"
test_that("Export IID for survival",
expect_equal(apply(resIID$survival.iid, 1:2, function(x){sqrt(sum(x^2))}),
             resSE$survival.se)
)

#### Cox model ####
iidTEST <-  iidCox(coxph.fit)

# }}}

# {{{ "cumsum(iid hazard) = iid cumhazard"
test_that("cumsum(iid hazard) = iid cumhazard",{
  Mcs <- t(apply(iidTEST$IFhazard[[1]], 1, cumsum))
  expect_equal(Mcs, iidTEST$IFcumhazard[[1]])
})

# }}}

# {{{ "iid equivalent parametrisation"
test_that("iid equivalent parametrisation",{
  expect_equal(iidTEST,
               iidCox(coxph.fit, newdata = d, tau.hazard = iidTEST$time[[1]])
  )
})

iidTEST2 <- iidCox(coxph.fit, tau.hazard = sort(c(iidTEST$time[[1]],c(0.1,1,5,8,12,25))))

# }}}
# {{{ "iid hazard = 0 at non event times"
test_that("iid hazard = 0 at non event times",{

  expect_equal(iidTEST2$IFhazard[[1]][,as.character(iidTEST$time[[1]])],
               iidTEST$IFhazard[[1]])

  expect_equal(iidTEST2$IFcumhazard[[1]][,as.character(iidTEST$time[[1]])],
               iidTEST$IFcumhazard[[1]])

  otherTimes0 <- setdiff(iidTEST2$time[[1]][iidTEST2$time[[1]]<max(coxph.fit$y[,1])],
                         iidTEST$time[[1]])
  expect_true(all(iidTEST2$IFhazard[[1]][,as.character(otherTimes0)] == 0))

  otherTimesNA <- iidTEST2$time[[1]][iidTEST2$time[[1]]>max(coxph.fit$y[,1])]
  expect_true(all(is.na(iidTEST2$IFhazard[[1]][,as.character(otherTimesNA)])))
})


iidTEST3 <- iidCox(coxph.fit, tau.hazard = iidTEST$time[[1]][2:3])

# }}}
# {{{ "iid hazard = 0 at non event times"
test_that("iid hazard = 0 at non event times",{
  
  expect_equal(iidTEST3$IFhazard[[1]],
               iidTEST$IFhazard[[1]][,as.character(iidTEST3$time[[1]])])
  
  expect_equal(iidTEST3$IFcumhazard[[1]],
               iidTEST$IFcumhazard[[1]][,as.character(iidTEST3$time[[1]])])
  
})

m.CSC <- CSC(Hist(time,event)~ X1+X2,data=d)

# }}}
# {{{ "iid ok when the last event is other cause"
test_that("iid ok when the last event is other cause",{
  d$status <- d$event==1
  IF1 <- iidCox(m.CSC$models$`Cause 1`, newdata = d[1:2,], tau.hazard = m.CSC$eventTimes)
  
  lastEvent <- m.CSC$models$`Cause 1`$y[tail(which(m.CSC$models$`Cause 1`$y[,2]==1),1),1]
  postLastEvent <- m.CSC$eventTimes[m.CSC$eventTimes>lastEvent]
  
  sapply(postLastEvent, function(t){
    expect_equal(IF1$IFcumhazard[[1]][,as.character(t)],
                 IF1$IFcumhazard[[1]][,as.character(lastEvent)])  
  })
  
})

# }}}

# {{{ "iid with newdata with an operator in status (here event==1 in Surv)"
test_that("iid with newdata with an operator in status (here event==1 in Surv)",{
    iidCox(coxph.fit, newdata = d)
})

# }}}





