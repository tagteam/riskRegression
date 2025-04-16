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
    CP <- c(0.1,0.3,0.8,0.89,0.9,0.91,7,28.0,29)
    get_SeSpPPVNPV <- function(cc){
        do.call(rbind,lapply(cc,function(c){
            suppressWarnings(res.SeSpPPVNPV.bili <- SeSpPPVNPV(cutpoint=c,T=pbc$time,delta=rep(1,nrow(pbc)),marker=pbc$bili,cause=1,weighting="marginal",times=1200,iid=TRUE))
            cbind(cutpoint = c, as.data.table(lapply(res.SeSpPPVNPV.bili[c(1,2,3,4)],"[",2)))
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
    x <- Score(list("Bilirubin" = pbc$bili),formula=Y~1,data=pbc,metrics="auc",breaks = sort(unique(pbc$bili)),cutpoints=CP)
    Score_result <- x$AUC$cutpoints[,.(cutpoint,TPR,FPR,PPV,NPV)]
    ROC_result <- x$AUC$ROC[,{
        pos <- 1+sindex(jump.times = risk,eval.times = CP)
        .(cutpoint = CP,TPR = c(0,TPR)[pos],FPR = c(0,FPR)[pos],PPV = c(1,PPV)[pos],NPV = c(0,NPV)[pos])
    }]

    ## Check that they give the same results
    expect_equal(as.numeric(x$AUC$cutpoints[,c("TPR","FPR","PPV","NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.bili[c(1,2,3,4)],"[",2)))
    expect_equal(as.numeric(x$AUC$cutpoints[,c("se.TPR","se.FPR","se.PPV","se.NPV")]),
                 as.numeric(sapply(res.SeSpPPVNPV.bili$inference[c(7,9,10,12)],"[",2)))
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
