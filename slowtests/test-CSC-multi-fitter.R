### test-CSC-multi-fitter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 25 2026 (10:14) 
## Version: 
## Last-Updated: mar 11 2026 (13:50) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
test_that("CSC supports fitter as list and fitter_arguments per model with ... override", {
    set.seed(9)
    d <- riskRegression::sampleData(200, outcome = "competing.risks")
    fit <- riskRegression::CSC(
                               formula = prodlim::Hist(time, event) ~ X1 + X2 + X8,
                               data = d,
                               cause = 1,
                               fitter = c("cph","coxph"),
                               fitter_arguments = list(list(method = "breslow"),list(ties = "breslow")),
                               ties = "efron"  
                           )
    # All models should be coxph fits because cph objects have both classes: cph and coxph
    expect_true(inherits(fit$models[[1]],"coxph"))
    expect_true(inherits(fit$models[[2]],"coxph"))
    # Check that ... overrode per-model ties (recorded in model$call)
    got_ties <- vapply(fit$models, function(m) {
        # ties stored in call for coxph
        as.character(m$call$ties)
    }, character(1))
    expect_true(all(got_ties == "efron"))
})

test_that("CSC errors if fitter list length mismatches number of models", {
    set.seed(9)
    d <- riskRegression::sampleData(120, outcome = "competing.risks")
    causes <- sort(unique(d$event[d$event > 0]))
    NC <- length(causes)
    expect_error(
        riskRegression::CSC(
                            prodlim::Hist(time, event) ~ X1 + X2,
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
    set.seed(9)
    d <- riskRegression::sampleData(200, outcome = "competing.risks")
    causes <- sort(unique(d$event[d$event > 0]))
    NC <- length(causes)
    obj <- riskRegression::CSC(
                               prodlim::Hist(time, event) ~ X1 + X2,
                               data = d,
                               cause = causes[1],
                               surv.type = "hazard",
                               fitter = rep(list("coxph"), NC),
                               fitter_arguments = list(list(ties = "efron"))
                           )

    got_ties <- vapply(obj$models, function(m) as.character(m$call$ties), character(1))
    expect_true(all(got_ties == "efron"))
})

test_that("CSC applies penalty.factor arguments correctly", {
    set.seed(9)
    d <- riskRegression::sampleData(200, outcome = "competing.risks")
    obj <- riskRegression::CSC(prodlim::Hist(time, event) ~ pen(X1,0) + unpenalized(X2) + pen(X8,8),data = d,cause = 1,fitter = c("glmnet","glmnet"),fitter_arguments = list(list(penalty.factor = c(X11 = 0,X21 = 1,X8 = 2))))
    expect_equal(obj$models[[1]]$penalty.factor,c(X11 = 0, X8 = 2, X21 = 1))
    expect_equal(obj$models[[2]]$penalty.factor,c(X11 = 0, X8 = 2, X21 = 1))
})



######################################################################
### test-CSC-multi-fitter.R ends here
