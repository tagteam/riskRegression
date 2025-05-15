library(testthat)
library(data.table)
library(survival)
library(timeROC)
data(Paquid)
data(pbc)
pbc<-pbc[!is.na(pbc$trt),] # select only randomised subjects
pbc$status<-as.numeric(pbc$status==2) # create event indicator: 1 for death, 0 for censored     

##-------------Without competing risks-------------------
test_that("pbc survival Sens,Spec,NPV,PPV",{
    # Se, Sp, PPV and NPV computation for serum bilirunbin at threshold c=0.9(mg/dl) 
    suppressWarnings(res.SeSpPPVNPV.bili <- SeSpPPVNPV(cutpoint=0.9,T=pbc$time,delta=pbc$status,marker=pbc$bili,cause=1,weighting="marginal",times=1200,iid=TRUE))
    x <- Score(list(pbc$bili),formula=Surv(time,status)~1,data=pbc,times=1200,metrics="auc",cutpoints=0.9)
    ## Check that they give the same results
    expect_equal(as.numeric(x$AUC$cutpoints[,c("TPR","FPR","PPV","NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.bili[c(1,2,3,4)],"[",2)))
    expect_equal(as.numeric(x$AUC$cutpoints[,c("se.TPR","se.FPR","se.PPV","se.NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.bili$inference[c(7,9,10,12)],"[",2)))
})
## Binary
test_that("pbc binary Sens,Spec,NPV,PPV",{
    pbc$Y <- 1*(pbc$time <= 1200)
    ## as.data.table(lapply(suppressWarnings(SeSpPPVNPV(cutpoint=0.8,T=pbc$time,delta=rep(1,nrow(pbc)),marker=pbc$bili,cause=1,weighting="marginal",times=1200,iid=TRUE))[c(1,2,3,4)],"[",2))
    ## as.data.table(lapply(suppressWarnings(SeSpPPVNPV(cutpoint=0.79,T=pbc$time,delta=rep(1,nrow(pbc)),marker=pbc$bili,cause=1,weighting="marginal",times=1200,iid=TRUE))[c(1,2,3,4)],"[",2))
    # the following cutpoints include values not in the data and outside the range
    CP <- c(0.1,0.3,0.79,0.8,0.81,0.89,0.9,0.91,7,28.0,29)
    get_SeSpPPVNPV <- function(cc){
        do.call(rbind,lapply(cc,function(c){
            suppressWarnings(x <- SeSpPPVNPV(cutpoint=c,T=pbc$time,delta=rep(1,nrow(pbc)),marker=pbc$bili,cause=1,weighting="marginal",times=1200,iid=TRUE))
            u = cbind(cutpoint = c,
                      as.data.table(lapply(x[c(1,2,3,4)],"[",2)),
                      as.data.table(lapply(x$inference[c(7,9,10,12)],"[",2)))
            setnames(u,c("cutpoint","TPR","FPR","PPV","NPV","se.TPR","se.FPR","se.PPV","se.NPV"))
            setcolorder(u,c("cutpoint","TPR","se.TPR","FPR","se.FPR","PPV","se.PPV","NPV","se.NPV"))
            u
        }))
    }
    manual_greater_or_equal <- function(pred,Y,cutpoints){
        do.call(rbind,lapply(cutpoints,function(c){
            data.table(cutpoint = c,
                       TPR = sum(pred >= c & Y == 1)/sum(Y == 1),
                       FPR = sum(pred >= c & Y == 0)/sum(Y == 0),
                       # if no one is predicted positive then the positive predictive rate is 1
                       PPV = if (sum(pred >= c) == 0) 1 else sum(pred >= c & Y == 1)/sum(pred >= c),
                       # if no one is predicted negative then the negative predictive rate is 1
                       NPV = if (sum(pred < c) == 0) 1 else sum(pred < c & Y == 0)/sum(pred < c))
        }))
    }
    # manual computation
    manual_result <- manual_greater_or_equal(pred = pbc$bili,Y = pbc$Y,cutpoints = CP)
    # timeROC
    timeROC_result <- get_SeSpPPVNPV(cc = CP)
    ## timeROC_result
    # Score
    x <- Score(list("Bilirubin" = pbc$bili),formula=Y~1,data=pbc,metrics="auc",breaks = sort(unique(pbc$bili)),cutpoints=CP)$AUC$cutpoints[,model := NULL][]
    ## rbind(x[cutpoint == 0.80][,model := NULL],timeROC_result[cutpoint == 0.80])
    ## rbind(x[cutpoint == 0.80][,model := NULL],timeROC_result[cutpoint == 0.80])
    ## Check that they give the same results (since timeROC uses < where score uses <= the results are only
    ## equal at the cutpoints which do not occur in the data.
    expect_equal(as.numeric(x[cutpoint == 0.80]), as.numeric(timeROC_result[cutpoint == 0.79]))
    expect_equal(as.numeric(x[cutpoint == 7]), as.numeric(timeROC_result[cutpoint == 7]))
})
test_that("Paquid competing risks Sens,Spec,NPV,PPV",{
    suppressWarnings(res.SeSpPPVNPV.DSST <- SeSpPPVNPV(cutpoint=22,
                                                       T=Paquid$time,
                                                       delta=Paquid$status,marker=Paquid$DSST,
                                                       cause=1,weighting="marginal",
                                                       times=5,iid=TRUE))
    x<-Score(list(Paquid$DSST),formula=Hist(time,status)~1,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
    ## Check that they give the same results
    expect_equal(as.numeric(x$AUC$cutpoints[,c("TPR","FPR","PPV","NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.DSST[c(1,3,4,6)],"[",2)))
    expect_equal(as.numeric(x$AUC$cutpoints[,c("se.TPR","se.FPR","se.PPV","se.NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.DSST$inference[c(7,9,10,12)],"[",2)))
    ## Test Cox censoring save memory
    x<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22)
    y<-Score(list(Paquid$DSST),formula=Hist(time,status)~DSST,data=Paquid,times=c(5),metrics="auc", cutpoints=22,censoring.save.memory = TRUE)
    expect_equal(x$AUC$cutpoints, y$AUC$cutpoints,tolerance=0.001)
    ## They are a different but this would usually not matter.
})
