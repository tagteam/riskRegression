### train.test for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
library(ranger)
context("riskRegression")

test_that("loob binary",{
    set.seed(8)
    learndat=sampleData(38,outcome="binary")
    f1 = glm(Y~X6,data=learndat,family = "binomial")
    f2 = glm(Y~X7+X8+X9,data=learndat,family = "binomial")
    ## leave-one-out bootstrap
    x <- Score(list(f1,f2),formula=Y~1,data=learndat,split.method="loob",B=100,se.fit=FALSE,seed = 5,verbose = -1,contrast = 0L,metrics = "auc")
    y <- Score(list(f1,f2),formula=Y~1,data=learndat,split.method="loob",B=10,se.fit=FALSE,seed = 5,verbose = -1,contrast = 0L,metrics = "auc")
})

test_that("loob survival",{
    set.seed(8)
    learndat=sampleData(38,outcome="survival")
    cox1a = coxph(Surv(time,event)~X6,data=learndat,x=TRUE,y=TRUE)
    cox2a = coxph(Surv(time,event)~X7+X8+X9,data=learndat,x=TRUE,y=TRUE)
    ## leave-one-out bootstrap
    x <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="loob",B=100,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
    y <- Score(list("COX1"=cox1a,"COX2"=cox2a),formula=Surv(time,event)~1,data=learndat,times=5,seed = 5,split.method="loob",B=10,se.fit=FALSE,progress.bar=NULL,metrics = "auc",verbose = -1)
})

#does give some warnings probably, nothing too serious
test.cv <- function(){
  type.cvs <- c("cv5","loob", "bootcv")
  B <- 100
  train.sample.size <- 400
  tau <- c(4,5)
  types <- list(list(data = "binary", model = "glm", response = "Y",left.extra = "family", right.extra = "binomial"),
                list(data = "survival", model = "coxph", response = "Surv(time,event)",left.extra = "x", right.extra = TRUE), 
                list(data = "competing.risks", model = "CSC", response = "Hist(time,event)"))
  for (cv in type.cvs){
      for (typ in types){
      cat(paste0("Testing ", typ$data, " data with ",cv, "  \n"))
      test_that(paste0("Testing ", typ$data, " data with ",cv),{
        set.seed(18)
        train.data <- sampleData(n=train.sample.size,outcome=typ$data)
        input.m <- list(data=train.data)
        if (typ$data != "competing.risks"){
          input.m[[typ$left.extra]] <- typ$right.extra
        }
        input.m1 <- append(list(formula = as.formula(paste0(typ$response,"~X1+X2+X7+X9"))), input.m) 
        input.m2 <- append(list(formula = as.formula(paste0(typ$response,"~X3+X5+X6"))), input.m)
        m1 <- do.call(typ$model, input.m1)
        m2 <- do.call(typ$model, input.m2)
        input.score <- list(object =list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2), formula = as.formula(paste0(typ$response,"~1")),data=train.data,conf.int=TRUE,progress.bar=NULL, split.method = cv, B=B)
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
}

test.cv()

test.cv2 <- function(){
  cv <- "cv5"
  B <- 1
  train.sample.size <- 400
  tau <- c(4,5)
  types <- list(list(data = "binary", model = "glm", response = "Y",left.extra = "family", right.extra = "binomial"),
                list(data = "survival", model = "coxph", response = "Surv(time,event)",left.extra = "x", right.extra = TRUE), 
                list(data = "competing.risks", model = "CSC", response = "Hist(time,event)"))
  for (typ in types){
    cat(paste0("Testing ", typ$data, " data with ",cv, "  \n"))
    test_that(paste0("Testing ", typ$data, " data with ",cv),{
      set.seed(18)
      train.data <- sampleData(n=train.sample.size,outcome=typ$data)
      input.m <- list(data=train.data)
      if (typ$data != "competing.risks"){
        input.m[[typ$left.extra]] <- typ$right.extra
      }
      input.m1 <- append(list(formula = as.formula(paste0(typ$response,"~X1+X2+X7+X9"))), input.m) 
      input.m2 <- append(list(formula = as.formula(paste0(typ$response,"~X3+X5+X6"))), input.m)
      m1 <- do.call(typ$model, input.m1)
      m2 <- do.call(typ$model, input.m2)
      input.score <- list(object =list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2), formula = as.formula(paste0(typ$response,"~1")),data=train.data,conf.int=TRUE,progress.bar=NULL, split.method = cv, B=B)
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

test.cv2()
