### train.test for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
context("riskRegression")

test.train.test <- function(){
  train.sample.size <- 200
  test.sample.size <- 500
  tau <- c(4,5)
  types <- list(list(data = "binary", model = "glm", response = "Y",left.extra = "family", right.extra = "binomial"),
                list(data = "survival", model = "coxph", response = "Surv(time,event)",left.extra = "x", right.extra = TRUE), 
                list(data = "competing.risks", model = "CSC", response = "Hist(time,event)"))
  for (typ in types){
    cat(paste0("Testing ", typ$data, " data with \n"))
    test_that(paste0("Testing ", typ$data, " data with "),{
      set.seed(18)
      train.data <- sampleData(n=train.sample.size,outcome=typ$data)
      test.data <- sampleData(n=test.sample.size,outcome=typ$data)
      input.m <- list(data=train.data)
      if (typ$data != "competing.risks"){
        input.m[[typ$left.extra]] <- typ$right.extra
      }
      input.m1 <- append(list(formula = as.formula(paste0(typ$response,"~X1+X2+X7+X9"))), input.m) 
      input.m2 <- append(list(formula = as.formula(paste0(typ$response,"~X3+X5+X6"))), input.m)
      m1 <- do.call(typ$model, input.m1)
      m2 <- do.call(typ$model, input.m2)
      input.score <- list(object =list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2), formula = as.formula(paste0(typ$response,"~1")),data=test.data,conf.int=TRUE,progress.bar=NULL)
      if (typ$data != "binary"){
        cat("with formula=Hist(time,event) ~ 1 and conservative = FALSE \n")
        input.score[["conservative"]] <- FALSE #check dependence on covariates
        input.score[["times"]] <- tau
        x1 <- do.call(Score, input.score)
        cat("with formula=Hist(time,event) ~ 1 and conservative = TRUE \n")
        input.score[["conservative"]] <- TRUE #check conservative when censoring does not depend on covariates
        x2 <- do.call(Score, input.score)
        cat("with formula=Hist(time,event) ~ X1+X2 and conservative = FALSE \n")
        input.score[["conservative"]] <- FALSE #check dependence on covariates
        input.score[["formula"]] <- as.formula(paste0(typ$response,"~X1+X2"))
        x3 <- do.call(Score, input.score)
        cat("with formula=Hist(time,event) ~ X1+X2 and conservative = TRUE \n")
        input.score[["conservative"]] <- TRUE #check dependence on covariates when conservative is true!
        x4 <- do.call(Score, input.score)
      }
      else {
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
}

test.train.test()
