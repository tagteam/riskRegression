### train.test for Score
library(testthat)
library(prodlim)
library(survival)
library(riskRegression)
library(data.table)
## library(ranger)
context("riskRegression")

## Testing binary data with cv5  
test_that("Testing binary data with cv5", {
    set.seed(18)
    train.data <- sampleData(n=50, outcome="binary")
    input.m <- list(data=train.data, family=binomial)
    input.m1 <- append(list(formula = Y ~ X1 + X2 + X7 + X9), input.m)
    input.m2 <- append(list(formula = Y ~ X3 + X5 + X6), input.m)
    m1 <- do.call(glm, input.m1)
    m2 <- do.call(glm, input.m2)
    input.score <- list(object = list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2), formula = Y ~ 1, data=train.data, conf.int=TRUE, progress.bar=NULL, split.method="cv5", B=100)
    x1 <- do.call(Score, input.score)
    expect_output(print(x1))
})

test_that("Testing survival data with cv5", {
    set.seed(18)
    train.data <- sampleData(n=70, outcome="survival")
    input.m <- list(data=train.data, x=TRUE)
    input.m1 <- append(list(formula = Surv(time, event) ~ X1 + X2 + X7 + X9), input.m)
    input.m2 <- append(list(formula = Surv(time, event) ~ X3 + X5 + X6), input.m)
    m1 <- do.call(coxph, input.m1)
    m2 <- do.call(coxph, input.m2)
    input.score <- list(object = list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2),
                        formula = Surv(time, event) ~ 1,
                        data=train.data,
                        conf.int=TRUE,
                        progress.bar=NULL,
                        split.method="cv5",
                        B=2)
    ## with formula=Hist(time,event) ~ 1 and conservative = FALSE 
    input.score$conservative <- FALSE
    input.score$times <- c(4,5)
    x1 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ 1 and conservative = TRUE 
    input.score$conservative <- TRUE
    x2 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ X1+X2 and conservative = FALSE 
    input.score$conservative <- FALSE
    input.score$formula <- Surv(time,event) ~ X1+X2
    x3 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ X1+X2 and conservative = TRUE 
    input.score$conservative <- TRUE
    x4 <- do.call(Score, input.score)
    expect_output(print(x1))
    expect_output(print(x2))
    expect_output(print(x3))
    expect_output(print(x4))
})

## Testing competing.risks data with cv5  
test_that("Testing competing.risks data with cv5", {
    set.seed(18)
    train.data <- sampleData(n=389, outcome="competing.risks")
    input.m <- list(data=train.data)
    input.m1 <- append(list(formula = Hist(time, event) ~ X1 + X7 + X9), input.m)
    input.m2 <- append(list(formula = Hist(time, event) ~ X1 + X6), input.m)
    m1 <- do.call(CSC, input.m1)
    m2 <- do.call(CSC, input.m2)
    input.score <- list(object = list("m(X1+X2+X7+X9)"=m1,"m(X3+X5+X6)"=m2),
                        formula = Hist(time, event) ~ 1,
                        data=train.data,
                        conf.int=TRUE,
                        progress.bar=NULL,
                        split.method="cv5",
                        B=2)
    ## with formula=Hist(time,event) ~ 1 and conservative = FALSE 
    input.score$conservative <- FALSE
    input.score$times <- c(4,5)
    x1 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ 1 and conservative = TRUE 
    input.score$conservative <- TRUE
    x2 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ X1+X2 and conservative = FALSE 
    input.score$conservative <- FALSE
    input.score$formula <- Hist(time,event) ~ X1+X2
    x3 <- do.call(Score, input.score)
    ## with formula=Hist(time,event) ~ X1+X2 and conservative = TRUE 
    input.score$conservative <- TRUE
    x4 <- do.call(Score, input.score)
    expect_output(print(x1))
    expect_output(print(x2))
    expect_output(print(x3))
    expect_output(print(x4))
})
