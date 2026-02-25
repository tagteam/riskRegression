### test-GLMnet-pen-rcs.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 25 2026 (10:14) 
## Version: 
## Last-Updated: feb 25 2026 (10:14) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
test_that("CSC supports fitter as list and fitter_arguments per model with ... override", {
  skip_if_not_installed("survival")
  skip_if_not_installed("prodlim")

  d <- riskRegression::sampleData(200, outcome = "competing.risks")

  # Ensure >1 cause present
  expect_true(length(unique(d$status[d$status > 0])) >= 2)

  # Determine causes and expected number of models for surv.type="hazard"
  causes <- sort(unique(d$status[d$status > 0]))
  NC <- length(causes)
  expect_true(NC >= 2)

  # Per-cause fitter list: both coxph but could differ; keep simple for CI stability
  fitters <- rep(list("coxph"), NC)

  # Different per-cause tie handling, and let ... override both to "breslow"
  fa <- vector("list", NC)
  fa[[1]] <- list(ties = "efron")
  if (NC >= 2) fa[[2]] <- list(ties = "exact")

  obj <- riskRegression::CSC(
    formula = prodlim::Hist(time, status) ~ X1 + X2,
    data = d,
    cause = causes[1],
    surv.type = "hazard",
    fitter = fitters,
    fitter_arguments = fa,
    ties = "breslow"  # should override per-cause ties
  )

  expect_s3_class(obj, "CauseSpecificCox")
  expect_type(obj$models, "list")
  expect_length(obj$models, NC)

  # All models should be coxph fits
  expect_true(all(vapply(obj$models, inherits, logical(1), "coxph")))

  # Check that ... overrode per-model ties (recorded in model$call)
  got_ties <- vapply(obj$models, function(m) {
    # ties stored in call for coxph
    as.character(m$call$ties)
  }, character(1))
  expect_true(all(got_ties == "breslow"))
})

test_that("CSC errors if fitter list length mismatches number of models", {
  skip_if_not_installed("prodlim")

  d <- riskRegression::sampleData(120, outcome = "competing.risks")
  causes <- sort(unique(d$status[d$status > 0]))
  NC <- length(causes)
  skip_if(NC < 2)

  expect_error(
    riskRegression::CSC(
      prodlim::Hist(time, status) ~ X1 + X2,
      data = d,
      cause = causes[1],
      surv.type = "hazard",
      fitter = list("coxph"),            # wrong length
      fitter_arguments = NULL
    ),
    "length must match|match the number of fitted models"
  )
})

test_that("CSC recycles fitter_arguments length 1 across models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("prodlim")

  d <- riskRegression::sampleData(200, outcome = "competing.risks")
  causes <- sort(unique(d$status[d$status > 0]))
  NC <- length(causes)
  skip_if(NC < 2)

  obj <- riskRegression::CSC(
    prodlim::Hist(time, status) ~ X1 + X2,
    data = d,
    cause = causes[1],
    surv.type = "hazard",
    fitter = rep(list("coxph"), NC),
    fitter_arguments = list(list(ties = "efron"))
  )

  got_ties <- vapply(obj$models, function(m) as.character(m$call$ties), character(1))
  expect_true(all(got_ties == "efron"))
})

######################################################################
### test-GLMnet-pen-rcs.R ends here
