library(riskRegression)
library(testthat)
library(rms)

test_that("cross-validated ROC curves", {
    d <- sampleData(173)
    v <- sampleData(227)
    f <- CSC(Hist(time,event)~X1+X8+X9+X10,data = d)
    x <- Score(list(f),data = v,formula = Hist(time,event)~1,plots = "roc")
    X <- Score(list(f),data = v,formula = Hist(time,event)~1,plots = "roc",split.method = "bootcv",B = 5)
    Y <- Score(list(f),data = v,formula = Hist(time,event)~1,plots = "roc",split.method = "cv5",B = 1)
    y <- Score(list(f),data = v,formula = Hist(time,event)~1,plots = "roc",split.method = "cv5",B = 2)
    Q <- Score(list(f),data = v,formula = Hist(time,event)~1,plots = "roc",split.method = "loob",B = 100)
    plotROC(x)
    plotROC(X,add = TRUE,col = 2)
    plotROC(Y,add = TRUE,col = 3)
    plotROC(y,add = TRUE,col = 4)
    plotROC(Q,add = TRUE,col = 5)
}
