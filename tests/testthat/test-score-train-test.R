### train.test for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
context("riskRegression")

# Helper function to run a score test for a given data type
run_score_test <- function(typ, test_description) {
  train.sample.size <- 200
  test.sample.size <- 500
  tau <- c(4, 5)

  test_that(test_description, {
    set.seed(18)
    train.data <- sampleData(n = train.sample.size, outcome = typ$data)
    test.data <- sampleData(n = test.sample.size, outcome = typ$data)
    input.m <- list(data = train.data)

    if (!is.null(typ$left.extra)) {
      input.m[[typ$left.extra]] <- typ$right.extra
    }

    input.m1 <- append(list(formula = as.formula(paste0(typ$response, "~X1+X2+X7+X9"))), input.m)
    input.m2 <- append(list(formula = as.formula(paste0(typ$response, "~X3+X5+X6"))), input.m)

    m1 <- do.call(typ$model, input.m1)
    m2 <- do.call(typ$model, input.m2)

    input.score <- list(
      object = list("m(X1+X2+X7+X9)" = m1, "m(X3+X5+X6)" = m2),                                                                                                   
      formula = as.formula(paste0(typ$response, "~1")),
      data = test.data,
      conf.int = TRUE,
      progress.bar = NULL
    )

    if (typ$data != "binary") {
      input.score[["conservative"]] <- FALSE
      input.score[["times"]] <- tau
      x1 <- do.call(Score, input.score)

      input.score[["conservative"]] <- TRUE
      x2 <- do.call(Score, input.score)

      input.score[["conservative"]] <- FALSE
      input.score[["formula"]] <- as.formula(paste0(typ$response, "~X1+X2"))
      x3 <- do.call(Score, input.score)

      input.score[["conservative"]] <- TRUE
      x4 <- do.call(Score, input.score)
    } else {
      x1 <- do.call(Score, input.score)
      x2 <- "NOT NEEDED HERE"
      x3 <- "NOT NEEDED HERE"
      x4 <- "NOT NEEDED HERE"
    }

    expect_output(print(x1))
    expect_output(print(x2))
    expect_output(print(x3))
    expect_output(print(x4))
  })
}

# Function to test binary data
test_binary_data <- function() {
  run_score_test(
    list(data = "binary", model = "glm", response = "Y", left.extra = "family", right.extra = "binomial"),
    "Binary outcome"
  )
}

# Function to test survival data
test_survival_data <- function() {
  run_score_test(
    list(data = "survival", model = "coxph", response = "Surv(time,event)", left.extra = "x", right.extra = TRUE),
    "Survival outcome"
  )
}

# Function to test competing risks data
test_competing_risks_data <- function() {
  run_score_test(
    list(data = "competing.risks", model = "CSC", response = "Hist(time,event)"),
    "Competing risks outcome"
  )
}

# Run the tests
test_binary_data()
test_survival_data()
test_competing_risks_data()