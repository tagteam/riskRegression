### test-predictCox-basehaz.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (10:38) 
## Version: 
## last-updated: Mar 26 2018 (07:57) 
##           By: Thomas Alexander Gerds
##     Update #: 22
#----------------------------------------------------------------------
## 
### Commentary: 
##
## Compare the output of riskRegression:::predictCox with survival:::basehaz.
## In particular check that they output the same results when centered = TRUE and centered = FALSE.
## Also perform additional tests on the other argument of predictCox.
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(testthat)
library(rms)
library(survival)
# {{{ 1- [predictCox] Check baseline hazard 
cat("Estimation of the baseline hazard \n")
data(Melanoma)
# {{{ "baseline hazard - match basehaz results"
test_that("baseline hazard - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + invasion + ici, data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
               basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumhazard, 
               basehaz(fit.cph)$hazard, tolerance = 1e-8)
  
  ## possible differences due to different fit - coef(fit.coxph)-coef(fit.cph)
  expect_equal(predictCox(fit.cph),
               predictCox(fit.coxph, centered = TRUE), 
               tolerance = 100*max(abs(coef(fit.coxph)-coef(fit.cph))))
})

#### strata

# }}}
# {{{ "baseline hazard (strata) - order of the results"
test_that("baseline hazard (strata) - order of the results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, keep.strata = TRUE)$strata, 
               basehaz(fit.coxph)$strata)
  expect_equal(predictCox(fit.cph, keep.strata = TRUE)$strata, 
               basehaz(fit.cph)$strata)
  # expect_equal(as.numeric(predictCox(fit.coxph, keep.strata = TRUE)$strata), 
  #              as.numeric(predictCox(fit.cph, keep.strata = TRUE)$strata))
})

# }}}
# {{{ "baseline hazard (strata) - match basehaz results"
test_that("baseline hazard (strata) - match basehaz results",{
  fit.coxph <- coxph(Surv(time,status == 1) ~ thick + strata(invasion) + strata(ici), data = Melanoma, y = TRUE, x = TRUE)
  fit.cph <- cph(Surv(time,status == 1) ~ thick + strat(invasion) + strat(ici), data = Melanoma, y = TRUE, x = TRUE)
  
  expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
               basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
  expect_equal(predictCox(fit.cph)$cumhazard, 
               basehaz(fit.cph)$hazard)
  
})

# }}}
# }}}
# {{{ 2- [predictCox] Check baseline hazard with time varying variables
## example from help(coxph)

test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 

# {{{ "baseline hazard (time varying) - match basehaz results"
test_that("baseline hazard (time varying) - match basehaz results",{
    fit.coxph <- coxph(Surv(start, stop, event) ~ x, data = test2, x = TRUE, y = TRUE)
    fit.cph <- cph(Surv(start, stop, event) ~ x, data = test2, x = TRUE, y = TRUE)

    expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
                 basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
                 basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.cph)$cumhazard, 
                 basehaz(fit.cph)$hazard, tolerance = 1e-8)
  
})

test2.strata <- rbind(cbind(as.data.table(test2),S = 1),
                      cbind(as.data.table(test2),S = 2))
test2.strata$randomS <- rbinom(NROW(test2.strata),size = 1, prob = 1/2)

# }}}
# {{{ "baseline hazard (strata, time varying) - match basehaz results"
test_that("baseline hazard (strata, time varying) - match basehaz results",{
    fit.coxph <- coxph(Surv(start, stop, event) ~ strata(S) + x, data = test2.strata, x = TRUE, y = TRUE)
    fit.cph <- cph(Surv(start, stop, event) ~ strat(S) + x, data = test2.strata, x = TRUE, y = TRUE)

    expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
                 basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
                 basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.cph)$cumhazard, 
                 basehaz(fit.cph)$hazard, tolerance = 1e-8)
    
    fit.coxph <- coxph(Surv(start, stop, event) ~ strata(randomS) + x, data = test2.strata, x = TRUE, y = TRUE)
    fit.cph <- cph(Surv(start, stop, event) ~ strat(randomS) + x, data = test2.strata, x = TRUE, y = TRUE)

    expect_equal(predictCox(fit.coxph, centered = FALSE)$cumhazard, 
                 basehaz(fit.coxph, centered = FALSE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.coxph, centered = TRUE)$cumhazard, 
                 basehaz(fit.coxph, centered = TRUE)$hazard, tolerance = 1e-8)
    expect_equal(predictCox(fit.cph)$cumhazard, 
                 basehaz(fit.cph)$hazard, tolerance = 1e-8)

})
# }}}
# }}}


#----------------------------------------------------------------------
### test-predictCox-basehaz.R ends here
