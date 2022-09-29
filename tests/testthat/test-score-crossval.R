### test crossval/bootcv/loob for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
context("riskRegression")

test.cross <- function(){
  sample.size <- 200
  sample.size.comp <- 400 #competing risks need larger data sets for convergence
  B <- 100
  B.cv5 <- 5 # should be smaller than the other cases; otherwise it is too slow!
  M <- 0.8 #to be certain of convergence in the case of competing risks, this should be set high!
  tau <- c(4,5)
  types <- list(list(data = "binary", model = "glm", response = "Y",left.extra = "family", right.extra = "binomial", sample.size = sample.size),
                list(data = "survival", model = "coxph", response = "Surv(time,event)",left.extra = "x", right.extra = TRUE, sample.size = sample.size), 
                list(data = "competing.risks", model = "CSC", response = "Hist(time,event)", sample.size = sample.size.comp))
  split.methods <- list(list(name = "cv5", B=B.cv5 ) ,list(name="loob", B=B), list(name= "bootcv", B=B))
  for (typ in types){
    for (split in split.methods){
      cat(paste0("Testing ", typ$data, " data with ", split$name, "\n"))
      test_that(paste0("Testing ", typ$data, " data with ", split$name),{
      set.seed(18)
      train.data <- sampleData(n=typ$sample.size,outcome=typ$data)
      input.m <- list(data=train.data)
      if (typ$data != "competing.risks"){
        input.m[[typ$left.extra]] <- typ$right.extra
      }
      input.m1 <- append(list(formula = as.formula(paste0(typ$response,"~X1+X2+X7+X9"))), input.m) 
      input.m2 <- append(list(formula = as.formula(paste0(typ$response,"~X3+X5+X6"))), input.m)
      m1 <- do.call(typ$model, input.m1)
      m2 <- do.call(typ$model, input.m2)
      input.score <- list(object =list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2), formula = as.formula(paste0(typ$response,"~1")),data=train.data,conf.int=TRUE,split.method=split$name,B=split$B,progress.bar=NULL)
      if (split$name!="cv5"){
        input.score[["M"]] <- M
      }
      if (typ$data != "binary"){
        cat("with formula=Hist(time,event) ~ 1 and conservative = FALSE \n")
        input.score[["conservative"]] <- FALSE #check dependence on covariates
        input.score[["times"]] <- tau
      }
      x1 <- do.call(Score, input.score)
      if (typ$data != "binary"){
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
}

test.cross()

