library(riskRegression)
library(rms)
library(survival)
library(testthat)
library(data.table)
# {{{ header with data  
set.seed(10)
n <- 5e1
dt <- sampleData(n,outcome="competing.risks")
dt[,status:=1*(event!=0)]
dt$time.ties <- round(dt$time)
dt$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
handler <- if (Sys.info()["sysname"] == "Windows") "foreach" else "mclapply"
verbose <- FALSE
# }}}
