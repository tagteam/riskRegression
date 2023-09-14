### test-plotRisk.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 15 2022 (16:04)
## Version:
## Last-Updated: Sep 16 2022 (12:56)
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
library(riskRegression)
library(testthat)
library(rms)
library(prodlim)
library(cmprsk)
library(survival)
library(data.table)
library(lava)
test_that("More than one competing risk", {
  set.seed(8)
  learndat <- sampleData(80, outcome = "competing.risk")
  testdat <- sampleData(140, outcome = "competing.risk")
  setkey(learndat, time)
  setkey(testdat, time)
  learndat[, event := as.character(event)]
  testdat[, event := as.character(event)]
  learndat[9:17, event := "cr2"]
  testdat[9:17, event := "cr2"]
  m1 <- FGR(Hist(time, event) ~ X2 + X7 + X9, data = learndat, cause = 1)
  m2 <- CSC(Hist(time, event) ~ X2 + X7 + X9, data = learndat, cause = 1)
  xcr <- Score(list("FGR" = m1, "CSC" = m2),
    formula = Hist(time, event) ~ 1,
    data = testdat, summary = "risks", null.model = 0L, times = c(3, 5)
  )
  plotRisk(xcr, times = 3)
  # check when no censored before time horizon
  testdat[time <= 3 & event == 0, event := 1]
  xcr <- Score(list("FGR" = m1, "CSC" = m2),
    formula = Hist(time, event) ~ 1,
    data = testdat, summary = "risks", null.model = 0L, times = c(3, 5)
  )
  plotRisk(xcr, times = 3)
  # check when no censored in all data
  testdat[event == 0, event := 1]
  xcr <- Score(list("FGR" = m1, "CSC" = m2),
    formula = Hist(time, event) ~ 1,
    data = testdat, summary = "risks", null.model = 0L, times = c(3, 5)
  )
  plotRisk(xcr, times = 3)
  # check when no event-free at horizon
  testdat[event == 0, event := 1]
  Score(list("FGR" = m1, "CSC" = m2),
    formula = Hist(time, event) ~ 1,
    data = testdat, summary = "risks", null.model = 0L, times = c(3, 8)
  )
  # all predicted risks of model m2 are NA
  expect_error(xcr = Score(list("FGR" = m1, "CSC" = m2),
    formula = Hist(time, event) ~ 1,
    data = testdat, summary = "risks", null.model = 0L, times = c(3, 8)
  ))
})


######################################################################
### test-plotRisk.R ends here
