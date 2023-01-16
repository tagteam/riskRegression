### train.test for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
library(ranger)
context("riskRegression")

# 10-fold cross-validation
set.seed(10)
d <- sampleData(2070,outcome = "binary")
val <- sampleData(100000,outcome = "binary")
d[,Y:=factor(Y,levels=c("0","1"),labels=c("0","1"))]
f1 <- glm(Y~X1+X2+X8+X9,data = d,family = "binomial")
f2 <- glm(Y~X3+X4+X5+X6+X10,data = d,family = "binomial")
f3 <- ranger(Y~X3+X4+X5+X6+X10,data = d)
# --------------------------------------------------------------
# in large validation set
# --------------------------------------------------------------
xval <- Score(list(f1,f2,f3),data = val,formula = Y~1)
## > xval$Brier$score
        ## model     Brier           se     lower     upper
## 1: Null model 0.2479651 0.0001420693 0.2476866 0.2482435
## 2:        glm 0.2134582 0.0005520726 0.2123762 0.2145403
## 3:      glm.1 0.1972720 0.0005895840 0.1961165 0.1984276
## 4:     ranger 0.2147243 0.0008084981 0.2131397 0.2163089
# --------------------------------------------------------------
# loob estimates are similar
# --------------------------------------------------------------
loob <- Score(list(f1,f2,f3),data = d,split.method = "loob",B = 200,formula = Y~1)
        ## model     Brier           se     lower     upper
## 1: Null model 0.2492589 0.0007510435 0.2477869 0.2507309
## 2:        glm 0.2084237 0.0037779047 0.2010191 0.2158282
## 3:      glm.1 0.1951589 0.0041166662 0.1870904 0.2032274
## 4:     ranger 0.2181784 0.0053378286 0.2077164 0.2286403
# --------------------------------------------------------------
# but cv-k gives completely wrong results
# --------------------------------------------------------------
x10 <- Score(list(f1,f2,f3),data = d,split.method = "cv10",formula = Y~1)
x5 <- Score(list(f1,f2,f3),data = d,split.method = "cv5",formula = Y~1)
## > x5$Brier$score
        ## model      Brier           se      lower      upper
## 1: Null model 0.06226608 0.0003806768 0.06151996 0.06301219
## 2:        glm 0.05189180 0.0018768193 0.04821330 0.05557030
## 3:      glm.1 0.04862532 0.0020616343 0.04458459 0.05266605
## 4:     ranger 0.05341315 0.0027774230 0.04796950 0.05885680

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
